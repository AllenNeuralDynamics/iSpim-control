"""Abstraction of Mesospim-Prime instrument."""

import numpy as np  # for live view
import cv2
from datetime import timedelta, datetime
import calendar
import tifffile
from mock import NonCallableMock as Mock
from threading import Thread, Event
from collections import deque
from pathlib import Path
from math import ceil
from time import perf_counter, sleep
from .devices.dcamapi4 import *
from .devices.dcam import Dcam, Dcamapi
from .devices.ni import WaveformHardware
from .devices.tiger_components import SamplePose, CameraPose, FilterWheel
from .mesospim_config import MesospimConfig, MAX_OPEN_LOOP_Z_DIST_UM
from .spim_base import Spim
from .tiff_transfer import TiffTransfer
from . import compute_waveforms
# Drivers in external packages
from tigerasi.tiger_controller import TigerController, UM_TO_STEPS
from tigerasi.sim_tiger_controller import TigerController as SimTiger

# Mapping from configuration name to class
LASER_MODEL_TO_OBJ = \
    {
       "Vortran": Vortran,
        "ObisLS": ObisLS,
        "Oxxius": Oxxius # NEED TO ADD OXXIUS DRIVER HERE
    }

# Constants
IMG_MIN = 0
IMG_MAX = 65535


class Mesospim(Spim):

    def __init__(self, config_filepath: str,
                 log_filename: str = 'debug.log',
                 console_output: bool = True,
                 color_console_output: bool = False,
                 console_output_level: str = 'info',
                 simulated: bool = False):
        """Read config file. Create Mesopim components according to config."""
        super().__init__(config_filepath, log_filename, console_output,
                         color_console_output, console_output_level, simulated)
        # Config
        self.cfg = MesospimConfig(config_filepath)

        # Thread handles and objects for syncing.
        self.image_capture_worker = None  # captures images during volumetric
                                          # image capture sequence.
        self.live_view_worker = None  # for displaying a window of the current
                                      # image.
        self.livestream_worker = None  # captures images during livestream
        self.livestream_enabled = Event()
        self.live_view_enabled = Event()
        self.image_in_hw_buffer = Event()
        self.img_deque = deque(maxlen=2)  # circular buffer

        self.autoexpose = False
        self.img_limits = [IMG_MIN, IMG_MAX]

        # Setup Hardware Components.
        if not (self.simulated or Dcamapi.init()):
            error_msg = "failed to initialize Dcamapi. Is the " \
                        "hardware turned-on/plugged-in? Is the dcamapi.dll " \
                        "(windows) or libdcamapi.so (linux) library " \
                        "installed? Is another program communicating with " \
                        "the camera?"
            self.log.error(error_msg)
            raise RuntimeError(error_msg)
        device_status_text = "Simulating" if self.simulated else "Initializing"
        self.log.info(f"{device_status_text} external devices.")
        self.cam = {} # create list of dcam objects
        for n in self.cfg.ncams: # loop over cameras
            self.cam.append(Dcam() if not self.simulated else Mock(Dcam))
        self.tigerbox = TigerController(**self.cfg.tiger_obj_params) if not \
            self.simulated else SimTiger(**self.cfg.tiger_obj_params)
        self.sample_pose = SamplePose(self.tigerbox)
        self.camera_pose = CameraPose(self.tigerbox)
        self.filter_wheel = FilterWheel(self.tigerbox,
                                        **self.cfg.filter_wheel_params)
        self.daq = WaveformHardware(**self.cfg.daq_obj_params) if not \
            self.simulated else Mock(WaveformHardware)
        self.lasers = {}  # populated in _setup_lasers.

        # Internal state attributes.
        self.active_laser = None  # Bookkeep which laser is configured.
        self.live_status = None  # Bookkeep if we are running in live mode.

        # Apply config-specific configurations to each component.
        self._setup_camera()
        self._setup_etl()
        self._setup_lasers()
        # Record the software state.
        self._log_git_hashes()

        # Extra Internal State attributes for the current image capture
        # sequence. These really only need to persist for logging purposes.
        self.image_index = 0  # current image to capture.
        self.total_tiles = 0  # tiles to be captured.
        self.stage_y_pos = None
        self.stage_z_pos = None

    def configure_daq(self, active_wavelenth: int, live: bool = False):
        """Setup DAQ to play waveforms according to wavelength to image with.

        Note: DAQ waveforms depend on the laser wavelength we are imaging with.

        :param active_wavelenth: the laser wavelength to configure the daq for.
            Only one laser wavelength may be turned on at a time.
        :param live: True if we want the waveforms to play on loop.
        """
        # TODO: cache this.
        self.log.debug(f"Computing waveforms for {active_wavelenth}[nm] laser.")
        _, voltages_t = compute_waveforms.generate_waveforms(self.cfg,
                                                             active_wavelenth)
        self.log.debug("Setting up daq.")
        period_time = self.cfg.get_daq_cycle_time(active_wavelenth)
        ao_names_to_channels = self.cfg.daq_ao_names_to_channels
        self.daq.configure(period_time, ao_names_to_channels, live)
        self.daq.assign_waveforms(voltages_t)
        # TODO: Consider plotting waveforms and saving to data collection
        #  folder at the beginning of an imaging run.

    def _setup_camera(self):
        """Setup camera for lightsheet scanning with specific config params."""
        self.log.debug("Setting up camera.")
        if self.simulated: # Can't do much better here in terms of simulating.
            self.cam.buf_getlastframedata.return_value = np.zeros((self.cfg.camx, self.cfg.camy), dtype='uint16')
            self.cam.wait_capevent_frameready.return_value = True
            return
        for n in self.cfg.ncams:
            self.cam[n].dev_open()  # open camera
            # Set sensor mode to lightsheet.
            self.cam[n].prop_setvalue(DCAM_IDPROP.SENSORMODE,
                                   DCAMPROP.SENSORMODE.PROGRESSIVE)
            # Set line interval.
            self.cam[n].prop_setvalue(DCAM_IDPROP.INTERNAL_LINEINTERVAL,
                                   self.cfg.row_interval)
            # Set any-pixel on-time.
            self.cam[n].prop_setvalue(DCAM_IDPROP.EXPOSURETIME,
                                   self.cfg.row_exposure_time)
            # Configure external trigger.
            self.cam[n].prop_setvalue(DCAM_IDPROP.TRIGGERSOURCE,
                                   DCAMPROP.TRIGGERSOURCE.EXTERNAL)
            # Set lightsheet scan direction.
            # match this to galvo ramp direction, can be "FORWARD" or "BACKWARD"
            readout_setting = DCAMPROP.READOUT_DIRECTION[self.cfg.scan_direction]
            self.cam[n].prop_setvalue(DCAM_IDPROP.READOUT_DIRECTION, readout_setting)

            # Poll the camera directly to get the settings we wrote. Log them.
            self.log.debug("SENSORMODE setting: "
                           f"{self.cam[n].prop_getvalue(DCAM_IDPROP.SENSORMODE)}")
            self.log.debug("INTERNAL_LINEINTERVAL setting: "
                           f"{self.cam[n].prop_getvalue(DCAM_IDPROP.INTERNAL_LINEINTERVAL):.3f}")
            self.log.debug("EXPOSURETIME setting: "
                           f"{self.cam[n].prop_getvalue(DCAM_IDPROP.EXPOSURETIME):.3f}")
            self.log.debug("Trigger Global Exposure setting: "
                           f"{str(self.cam[n].prop_getvalue(DCAM_IDPROP.TRIGGER_GLOBALEXPOSURE))}")
            self.log.debug("Internal Line Speed setting: "
                           f"{str(self.cam[n].prop_getvalue(DCAM_IDPROP.INTERNALLINESPEED))}")

    def _setup_lasers(self):
        """Setup lasers that will be used for imaging. Warm them up, etc."""
        for wavelength_str, specs in self.cfg.laser_specs.items():
            self.log.debug(f"Setting up {specs['color']} laser.")
            cls = LASER_MODEL_TO_OBJ[specs['driver']]
            # Mock Lasers. We can't wrap them in Mock(cls(**specs['kwds'])
            laser_obj = cls(**specs['kwds']) if not self.simulated else Mock(cls)
            # Put all lasers in analog mode.
            # TODO: can this be abstracted, not case-specific?
            if cls == ObisLS:
                self.log.debug("Applying analog mode settings to ObisLS laser.")
                laser_obj.set_modulation_mode(LSModulationType.ANALOG)
                laser_obj.set_analog_input_impedance(
                    AnalogInputImpedanceType.TWO_THOUSAND_OHM)
            self.lasers[int(wavelength_str)] = laser_obj

    def _log_git_hashes(self):
        """Log the git hashes of this project and all packages."""
        pass

    def run_from_config(self):
        """Collect volumetric image according to the config file parameters."""
        if self.livestream_enabled.is_set():
            self.stop_livestream()
        self.capture_tiled_image_stack(self.cfg.volume_x_um,
                                       self.cfg.volume_y_um,
                                       self.cfg.volume_z_um,
                                       self.cfg.imaging_wavelengths,
                                       self.cfg.tile_overlap_y_percent,
                                       self.cfg.tile_overlap_z_percent,
                                       self.cfg.tile_prefix,
                                       self.cfg.local_storage_dir,
                                       self.img_storage_dir)

    def capture_tiled_image_stack(self, volume_x_um: float, volume_y_um: float,
                                  volume_z_um: float,
                                  wavelengths: list,
                                  tile_overlap_y_percent: float = 15.0,
                                  tile_overlap_z_percent: float = 15.0,
                                  tile_prefix: str = "tile",
                                  local_storage_dir: Path = Path("."),
                                  img_storage_dir: Path = None):
        """ Capture an array of tiles starting from the current position.
        Note: returns to the starting position even if errors occur
              (except stage-related errors).
        """
        # Iterate through the volume through z, then x, then y.
        # Play waveforms for moving the laser and sensor shutter.
        # Capture the fully-formed image after waveform has finished playing.
        # Move the stage by the specified amount between image capture events.

        # micrometers per grid step. At 0 tile overlap, this is just the
        # sensor's field of view.
        y_grid_step_um = \
            (1 - tile_overlap_y_percent/100.0) * self.cfg.tile_size_um
        z_grid_step_um = \
            (1 - tile_overlap_z_percent/100.0) * self.cfg.tile_size_um

        # Compute step count.
        # Always round up so that we cover the desired imaging region.
        xsteps = ceil((volume_z_um - self.cfg.x_voxel_size_um)
                      / self.cfg.x_voxel_size_um)
        ysteps = ceil((volume_y_um - self.cfg.tile_size_um)
                      / y_grid_step_um)
        zsteps = ceil((volume_z_um - self.cfg.tile_size_um)
                      / z_grid_step_um)

        # TODO: check if we will violate storage directory limits with our
        #   run. (Check local and external storage.)

        # Log relevant info about this imaging run.
        self.total_tiles = (1+xsteps)*(1+ysteps)*(1+zsteps)*len(wavelengths)
        self.log.info(f"Total tiles: {self.total_tiles}.")
        run_time_days = self.total_tiles/(self.cfg.tiles_per_second*3600.*24.)
        completion_date = datetime.now() + timedelta(days=run_time_days)
        date_str = completion_date.strftime("%d %b, %Y at %H:%M %p")
        weekday = calendar.day_name[completion_date.weekday()]
        self.log.info(f"Time Esimate: {run_time_days:.2f} days. Imaging run "
                      f"should finish after {weekday}, {date_str}.")
        self.log.info(f"Desired dimensions: {volume_x_um:.1f}[um] x "
                      f"{volume_y_um:.1f}[um] x {volume_z_um:.1f}[um]")
        self.actual_vol_x_um = self.cfg.x_voxel_size_um*(1+xsteps)
        self.actual_vol_y_um = self.cfg.tile_size_um + xsteps*y_grid_step_um
        self.actual_vol_z_um = self.cfg.tile_size_um + ysteps*z_grid_step_um
        self.log.info(f"Actual dimensions: {actual_vol_x_um:.1f}[um] x "
                      f"{actual_vol_y_um:.1f}[um] x {actual_vol_z_um:.1f}[um]")
        self.log.info(f"X grid step: {self.cfg.x_voxel_size_um} [um]")
        self.log.info(f"Y grid step: {y_grid_step_um} [um]")
        self.log.info(f"Z grid step: {z_grid_step_um} [um]")

        # Set the sample starting location as the origin.
        self.sample_pose.home_in_place()
        # Disable backlash compensation on X to improve x-axis scanning.
        # Disabling X compensation is ok since we scan in one direction.
        self.tigerbox.set_axis_backlash(x=0.0)
        # Apply a lead-in-move to take out backlash.
        x_backup_pos = -UM_TO_STEPS*self.cfg.stage_backlash_reset_dist_um
        self.log.debug("Applying extra move to take out backlash.")
        self.sample_pose.move_absolute(x=round(x_backup_pos), wait=True)
        self.sample_pose.move_absolute(x=0, wait=True)
        max_intensity_worker = None
        transfer_process = None
        self.stage_y_pos, self.stage_z_pos = (0, 0)
        self.image_index = 1
        # Setup double buffer to capture images while reading out prior image.
        self.cam.buf_alloc(2)
        self.cam.cap_start()
        try:
            for j in range(0, zsteps + 1):
                self.stage_x_pos = 0
                self.sample_pose.move_absolute(y=round(self.stage_z_pos),
                                               wait=True)
                for i in range(0, ysteps + 1):
                    self.sample_pose.move_absolute(x=round(self.stage_y_pos),
                                                   wait=True)
                    # Log temperature for this XY start location.
                    # self.log.info(f"etl_temp: {self.etl.temp_reading()}") change this to log through Tiger
                    # Collect a zstack for every specified laser/filter combo.
                    for wavelength in wavelengths:
                        # Setup capture of next Z stack.
                        # Open filters, enable active laser; disable the rest.
                        self.setup_imaging_for_laser(wavelength)
                        filename = Path(f"{tile_prefix}_{i}_{j}_{wavelength}.tiff")
                        filepath_src = local_storage_dir / filename
                        image_count = zsteps + 1
                        self._collect_stacked_tiff(image_count, wavelength,
                                                   filepath_src)
                        # Start transferring tiff file to its destination.
                        # Note: Image transfer is faster than image capture.
                        #   but we still wait for prior process to finish.
                        if transfer_process is not None:
                            self.log.info("Waiting for tiff transfer process "
                                          "to complete.")
                            # TODO: log how often this happens. Email someone
                            #  for help if it's happening a lot.
                            transfer_process.join()
                        self.log.info(f"Starting transfer process for {filename}.")
                        if ext_storage_dir is not None:
                            filepath_dest = ext_storage_dir / filename
                            transfer_process = TiffTransfer(filepath_src,
                                                            filepath_dest)
                            transfer_process.start()
                    self.stage_y_pos += y_grid_step_um * UM_TO_STEPS
                self.stage_z_pos += z_grid_step_um * UM_TO_STEPS
        finally:
            # Wait for waveform playback to finish so we leave signals in their
            # ending voltage states.
            self.daq.wait_until_done()  # use default timeout of 1[s].
            self.log.info("Stopping camera(s).")
            for n in self.cfg.ncams:
                self.cam[n].cap_stop()
                self.cam[n].buf_release()
            self.cam.buf_release()
            if transfer_process is not None:
                self.log.debug("joining zstack transfer process.")
                transfer_process.join()
            self.log.info("Returning to start position.")
            # Discard these queued replies before issuing move+wait.
            sleep(0.05)  # Hack to avoid race condition from any prior replies.
            self.tigerbox.clear_incoming_message_queue()
            self.sample_pose.move_absolute(x=0, y=0, z=0, wait=True)

    def _collect_stacked_tiff(self, image_count, wavelength: int,
                              filepath: Path):
        """helper fn to collect a stack of images and save it to file.

        This fn will close the file and save any images collected so far if
        aborted early.

        Ths fn assumes we are starting from x=0 and will leave the machine
        back at x=0. (A better way to do this would be to push/pop x location.)
        """
        # dispim scheme uses constant velocity stage scanning
        stage_x_pos = 0
        stage_x_backup_pos = -UM_TO_STEPS*self.cfg.stage_backlash_reset_dist_um

        images_captured = 0
        # setup the stage scanning for 1 line at stage_y_pos of length x_vol
        # need to look into how to disable auto-retrace.
        # either set z=0 in scanv, or configure f arg in scan command
        self.sample_pose.scanr(x=0, y=self.actual_vol_x_um/1000.0) # y is in units of mm
        self.sample_pose.scanv(x=self.stage_y_pos/1000.0, y=self.stage_y_pos/1000.0, z=1)
        self.sample_pose.ttl(x=1) # activate ttl interupt for encoder pulse output

        try:
            for k in range(0, image_count):
                loop_start = perf_counter()
                extras = {'x': self.stage_x_pos, 'y': self.stage_y_pos,
                          'z': stage_z_pos, 'image_index': self.image_index,
                          'monotonic_timestamp': loop_start}
                msg = "Snapping image " \
                      f"{self.image_index:12}/{self.total_tiles}. " \
                      f"sample position: ({self.stage_x_pos:.3f}, " \
                      f"{self.stage_y_pos:.3f}, {stage_z_pos:.3f}) [steps]"
                self.log.debug(msg, extra=extras)
                # Print this to the screen on the same line.
                if self.console_output:
                    print("\r" + msg, end=" ", flush=True)
                # Trigger waveforms to move laser and expose image.
                self.daq.start()
        finally:
            print()
            # If zwait was false, then we ignored many tigerbox replies.
            # Discard these queued replies before issuing move+wait.
            sleep(0.05)  # Hack to avoid race condition from any prior replies.
            self.tigerbox.clear_incoming_message_queue()
            # Reset Z stage position to where we started.
            # Apply a lead-in-move to take out backlash.
            self.log.debug("Applying extra move to take out backlash.")
            self.sample_pose.move_absolute(z=round(stage_z_backup_pos), wait=True)
            stage_z_pos = 0
            self.sample_pose.move_absolute(z=round(stage_z_pos), wait=True)
            # If we aborted early, save whatever pictures we took.
            if image_count != images_captured:
                # force tiff_capturer to timeout by setting image flag.
                self.image_in_hw_buffer.set()
            # Join the image capture thread.
            self.log.debug("Joining zstack capture thread.")
            tiff_capturer.join()

    # def _image_capture_worker(self, filepath: Path, image_count: int,
    #                           wavelength: int):
    #     images_captured = 0
    #     image_wait_time = round(5*self.cfg.get_daq_cycle_time(wavelength)*1e3)
    #     self.log.info(f"Creating tiff file: {filepath}")
    #     # Note: creating this file can take non-negligible time.
    #     with tifffile.TiffWriter(filepath,  bigtiff=True) as tif:
    #         while images_captured < image_count:
    #             capture_start = perf_counter()
    #             # Check if image readout is ready with reasonable timeout.
    #             self.image_in_hw_buffer.wait()
    #             if not self.cam.wait_capevent_frameready(image_wait_time):
    #                 raise RuntimeError(f"Error. Camera failed to capture image "
    #                                    f"{images_captured + 1} with error: "
    #                                    f"{str(self.cam.lasterr())}")
    #             self.img_deque.append(self.cam.buf_getlastframedata())
    #             tif.write(self.img_deque[-1], contiguous=True)
    #             self.image_in_hw_buffer.clear()
    #             images_captured += 1
    #             self.log.debug("Acquired zstack image: "
    #                            f"{images_captured}/{image_count}. "
    #                            f"Read took: {perf_counter() - capture_start:.3f}[s].")

    def setup_imaging_for_laser(self, wavelength: int, live: bool = False):
        """Configure system to image with the desired laser wavelength.

        Note: this fn changes the filter wheel too.
        """
        # Bail early if this laser is already setup in the previously set mode.
        if self.active_laser == wavelength and self.live_status == live:
            self.log.debug("Skipping daq setup. Laser already provisioned.")
            return
        self.live_status = live
        live_status_msg = " in live mode" if live else ""
        self.log.info(f"Configuring {wavelength}[nm] laser{live_status_msg}.")
        if self.active_laser is not None:
            self.lasers[self.active_laser].disable()
        fw_index = self.cfg.laser_specs[str(wavelength)]['filter_index']
        self.filter_wheel.set_index(fw_index)
        camera_focus_pos = self.cfg.laser_specs[str(wavelength)]['camera_axis_position']
        self.camera_pose.move_absolute(camera_focus_pos, wait=True)
        # Reprovision the DAQ.
        self.configure_daq(wavelength, live)
        self.active_laser = wavelength
        self.lasers[self.active_laser].enable()

    def get_latest_img(self):
        """returns the latest image as a 2d numpy array. Useful for UIs."""
        return self.img_deque[0]

    # def move_sample_absolute(self, x: int = None, y: int = None, z: int = None):
    #     """Convenience function for moving the sample from a UI."""
    #     self.sample_pose.move_absolute(x=x, y=y, z=z, wait=True)

    def move_sample_relative(self, x: int = None, y: int = None, z: int = None):
        """Convenience func for moving the sample from a UI (units: steps)."""
        self.sample_pose.move_relative(x=x, y=y, z=z, wait=True)

    def move_camera_relative(self, m: int):
        """Convenience func for moving the camera from a UI (units: steps)."""
        self.camera_pose.move_relative(m=m, wait=True)

    def move_camera_absolute(self, m: int):
        """Convenience func for moving the camera from a UI (units: steps)."""
        self.camera_pose.move_absolute(m, wait=True)

    def get_sample_position(self):
        return self.sample_pose.get_position()

    def get_camera_position(self):
        return self.camera_pose.get_position()

    def start_livestream(self, wavelength: int):
        """Repeatedly play the daq waveforms and buffer incoming images."""
        if wavelength not in self.cfg.laser_wavelengths:
            self.log.error(f"Aborting. {wavelength}[nm] laser is not a valid "
                           "laser.")
            return
        # Bail early if it's started.
        if self.livestream_enabled.is_set():
            self.log.warning("Not starting. Livestream is already running.")
            return
        self.log.debug("Starting livestream.")
        self.log.warning(f"Turning on the {wavelength}[nm] laser.")
        self.setup_imaging_for_laser(wavelength, live=True)
        self.daq.start()
        self.livestream_enabled.set()
        self.livestream_worker = Thread(target=self._livestream_worker,
                                        args=(wavelength,), daemon=True)
        self.livestream_worker.start()
        # Launch thread for picking up camera images.

    def stop_livestream(self, wait: bool = False):
        # Bail early if it's already stopped.
        if not self.livestream_enabled.is_set():
            self.log.warning("Not starting. Livestream is already stopped.")
            return
        wait_cond = "" if wait else "not "
        self.log.debug(f"Disabling livestream and {wait_cond}waiting.")
        self.livestream_enabled.clear()

        if wait:
            self.livestream_worker.join()
        self.daq.stop()
        self.daq.close()
        self.live_status = False  # TODO: can we get rid of this if we're always stopping the livestream?
        self.active_laser = None

    def _livestream_worker(self, wavelength):
        """Pulls images from the camera and puts them into the ring buffer."""
        image_wait_time = round(5*self.cfg.get_daq_cycle_time(wavelength)*1e3)
        self.cam.buf_alloc(2)
        self.cam.cap_start()
        while self.livestream_enabled.is_set():
            if not self.cam.wait_capevent_frameready(image_wait_time):
                self.log.error(f"Camera failed to capture image with error"
                               f"{str(self.cam.lasterr())}")
            if self.simulated:
                sleep(1./16)
            self.img_deque.append(self.cam.buf_getlastframedata())
        self.cam.cap_stop()
        self.cam.buf_release()


    def enable_display(self):
        """Enable live image viewing in a separate thread."""
        # Bail early if it's already enabled.
        if self.live_view_enabled.is_set():
            self.log.warning("Not enabling. Display is already enabled.")
            return
        self.log.debug("Enabling live view.")
        self.live_view_enabled.set()
        self.live_view_worker = Thread(target=self._live_view_worker,
                                       daemon=True)
        self.live_view_worker.start()

    def disable_display(self, wait: bool = False):
        if not self.live_view_enabled.is_set():
            self.log.warning("Not disabling. Display is already disabled.")
            return
        wait_cond = "" if wait else "not "
        self.log.debug(f"Disabling live view and {wait_cond}waiting.")
        self.live_view_enabled.clear()
        if wait:
            self.live_view_worker.join()

    def set_autoexpose(self, state: bool):
        """Set to True/False to Enable/Disable builtin display autoexposure."""
        self.autoexpose = state
        self.log.info(f"autoexpose is: {self.autoexpose} of type {type(self.autoexpose)}")

    def get_autoexposure(self):
        return self.img_limits

    def set_exposure_limits(self, min_limit: int = None,
                            max_limit: int = None):
        # Check extremes if defined.
        if min_limit is not None:
            if min_limit < IMG_MIN or min_limit > IMG_MAX:
                self.log.error("Invalid exposure minimum setting.")
                min_limit = self.img_limits[0]
            if max_limit < IMG_MIN or max_limit > IMG_MAX:
                self.log.error("Invalid exposure maximum setting")
                max_limit = self.img_limits[1]
        self.img_limits = [min_limit if min_limit is not None
                           else self.img_limits[0],
                           max_limit if max_limit is not None
                           else self.img_limits[1]]

    def _live_view_worker(self):
        """Threadworker for displaying a window to view the latest images."""
        window_name = "Live View"
        cv2.namedWindow(window_name, cv2.WINDOW_AUTOSIZE)
        frames_grabbed = 0
        while self.live_view_enabled.is_set():
            if len(self.img_deque) < 1:
                continue
            data = self.get_latest_img()
            # Shrink the data for the viewer and rotate the picture.
            data = np.rot90(cv2.resize(data, (640, 640)), -1)
            # Apply autoexposure if set.
            if self.autoexpose:
                self.img_limits[0] = IMG_MIN
                self.img_limits[1] = np.amax(data)
            # Apply current exposure settings.
            # Point Slope Formula.
            m = int(IMG_MAX/(self.img_limits[1] - self.img_limits[0]))
            data = np.round(np.clip(m*data - m*self.img_limits[0], IMG_MIN, IMG_MAX))
            cv2.imshow(window_name, data)
            key = cv2.waitKey(1)  # Display for at least one ms
        cv2.destroyAllWindows()

    def enable_live_mips(self):
        raise NotImplementedError

    def disable_live_mips(self):
        raise NotImplementedError

    def _live_mips_worker(self):
        """Threadworker for live computation of zstack max intensity projection."""
        # TODO: this might need to happen in a separate process with queued images
        #       so as not to slow down the main thread.
        # TODO: error if we cannot compute max intensity projection
        # fast enough to do it live. Bail but don't kill the run.
        raise NotImplementedError

    def reload_config(self):
        """Reload the toml file."""
        self.cfg.reload()

    def close(self):
        """Safely close all open hardware connections."""
        if self.livestream_enabled.is_set():
            self.stop_livestream()
        if self.live_view_enabled.is_set():
            self.disable_display()
        self.log.info("Closing ETL.")
        self.etl.close(soft_close=True)
        self.log.info("Closing daq.")
        self.daq.close()
        self.active_laser = None
        self.log.info("Powering down lasers.")
        for wavelength, laser in self.lasers.items():
            self.log.info(f"Powering down {wavelength}[nm] laser.")
            laser.disable()
        self.log.info("Closing camera.")
        self.cam.dev_close()
        self.log.info("De-initializing Dcam API.")
        # Dcamapi.uninit can be called repeatedly, but we must call it at least
        # once at the end.
        Dcamapi.uninit()
        super().close()