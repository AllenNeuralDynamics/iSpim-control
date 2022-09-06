#!/usr/bin/env python3
"""Abstraction of dispim instrument."""

import logging
import sys

from coloredlogs import ColoredFormatter
from datetime import timedelta, datetime
import calendar
import tifffile
from mock import NonCallableMock as Mock
from threading import Thread, Event
import os, shutil
import numpy as np
from pathlib import Path
from math import ceil
from time import perf_counter, sleep
from datetime import date
from .devices.dcamapi4 import *
from .devices.dcam import Dcam, Dcamapi
from .devices.ni import WaveformHardware
from .devices.tiger_components import SamplePose, CameraPose, FilterWheel
from .mesospim_config import MesospimConfig, MAX_OPEN_LOOP_Z_DIST_UM
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
        "Oxxius": Oxxisu # NEED TO OXXIUS ADD DRIVER HERE
    }


class Mesospim:

    def __init__(self, config_filepath: str,
                 log_filename: str = 'debug.log',
                 console_output: bool = True,
                 color_console_output: bool = False,
                 console_output_level: str = 'info',
                 simulated: bool = False):
        """Read config file. Create Mesopim components according to config."""
        # If simulated, no physical connections to hardware should be required.
        # Simulation behavior should be pushed as far down the object
        # hierarchy as possible.
        self.simulated = simulated
        # Setup logging.
        # Save console output to print/not-print imaging progress.
        self.console_output = console_output
        # We want the name of the package here since logger hierarchy
        # depends on module structure.
        self.log = logging.getLogger("dispim")
        # logger level must be set to the lowest level of any handler.
        self.log.setLevel(logging.DEBUG)
        # Create log handlers to dispatch:
        # - DEBUG level and above to write to a file called debug.log.
        # - User-specified level and above to print to console if specified.
        fmt = '%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s'
        fmt = "[SIM] " + fmt if self.simulated else fmt
        datefmt = '%Y-%m-%d,%H:%M:%S'
        self.log_format = logging.Formatter(fmt=fmt, datefmt=datefmt)
        self.log_handlers = []
        debug_filepath = Path(log_filename)
        self.log_handlers.append(logging.FileHandler(debug_filepath))
        self.log_handlers[-1].setLevel(logging.DEBUG)
        self.log_handlers[-1].setFormatter(self.log_format)
        if self.console_output:
            self.log_handlers.append(logging.StreamHandler(sys.stdout))
            self.log_handlers[-1].setLevel(console_output_level)
            if color_console_output:
                colored_formatter = ColoredFormatter(fmt=fmt, datefmt=datefmt)
                self.log_handlers[-1].setFormatter(colored_formatter)
            else:
                self.log_handlers[-1].setFormatter(self.log_format)
        for handler in self.log_handlers:
            self.log.addHandler(handler)

        # Config
        self.cfg = MesospimConfig(config_filepath)

        # Thread handles and objects for syncing.
        self.image_capture_worker = None
        self.live_view_worker = None
        self.live_view_enabled = Event()
        self.image_in_hw_buffer = Event()

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
        self.active_laser = None  # Bookkeeping which laser is configured.

        # Apply config-specific configurations to each component.
        self._setup_camera()
        self._setup_etl()
        self._setup_lasers()
        # Record the software state.
        self._log_git_hashes()

        # Extra attributes for the current volumetric image capture sequence.
        # These really only need to persist for logging purposes.
        self.image_index = 0  # current image to capture.
        self.total_tiles = 0  # tiles to be captured.
        self.stage_x_pos = None
        self.stage_y_pos = None

    def configure_daq(self, active_wavelenth: int):
        """Setup DAQ to play waveforms according to wavelength to image with.

        Note: DAQ waveforms depend on the laser wavelength we are imaging with.

        :param active_wavelenth: the laser wavelength to configure the daq for.
            Only one laser wavelength may be turned on at a time.
        """
        # TODO: cache this.
        self.log.debug(f"Computing waveforms for {active_wavelenth} [nm] laser.")
        _, voltages_t = compute_waveforms.generate_waveforms(self.cfg,
                                                             active_wavelenth)
        self.log.debug("Setting up daq.")
        period_time = self.cfg.get_daq_cycle_time(active_wavelenth)
        ao_names_to_channels = self.cfg.daq_ao_names_to_channels
        self.daq.configure(period_time, ao_names_to_channels)
        self.daq.assign_waveforms(voltages_t)
        # TODO: Consider plotting waveforms and saving to data collection
        #  folder at the beginning of an imaging run.

    def _setup_camera(self):
        """Setup camera for lightsheet scanning with specific config params."""
        self.log.debug("Setting up camera.")
        if self.simulated: # Can't do much better here in terms of simulating.
            self.cam.buf_getlastframedata.return_value = np.zeros((self.cfg.camx, self.cfg.camy), dtype='uint16')
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

    def livestream(self):
        """repeatedly play the daq waveform. Blocks.

        Can be called in conjunction with liveview to see the image.
        """
        self.log.debug("Starting livestream.")
        print("GUI settings for debugging:")
        print(f"  Line Interval: {self.row_interval*1e6:.3f} [us]")
        print(f"  Exposure Time (computed): {self.row_exposure_time*1e3:.3f} [ms]")
        print()
        self.enable_live_view()
        try:
            while True:
                self.daq.start()
                # TODO: possibly add a sleep here.
                self.daq.wait_until_done()  # use 1 [sec] default timeout.
                self.daq.stop()
        finally:
            self.log.debug("Stopping livestream.")
            self.disable_live_view()
            self.daq.wait_until_done()  # use 1 [sec] default timeout.
            self.daq.stop()

    def collect_volumetric_image(self, live_view: bool = False,
                                 compute_max_intensity_proj: bool = False,
                                 overwrite: bool = False):
        """collect volumetric image according to the config file parameters."""
        # Create output folder and folder for storing images.
        output_folder = \
            self.cfg.ext_storage_dir / Path(self.cfg.subject_id + "-ID_" +
                                            date.today().strftime("%Y_%m_%d"))
        print(f"external storage dir is: {self.cfg.ext_storage_dir.resolve()}")
        if output_folder.exists() and not overwrite:
            self.log.error(f"Output folder {output_folder.absolute()} exists, "
                           "This function must be rerun with overwrite=True.")
            raise
        img_folder = output_folder / Path("micr/")
        # Create a log file specific to this job.
        # TODO: enable/disable max proj computation in a separate process.
        imaging_log_filename = "imaging_log.log"  # This should be a cfg param.
        imaging_log_filepath = Path(imaging_log_filename)
        imaging_log_handler = logging.FileHandler(imaging_log_filepath)
        imaging_log_handler.setLevel(logging.DEBUG)
        imaging_log_handler.setFormatter(self.log_format)
        self.log.addHandler(imaging_log_handler)
        self._log_git_hashes()
        self.log.info(f"Creating datset folder in: {output_folder.absolute()}")
        try:
            print(img_folder)
            img_folder.mkdir(parents=True, exist_ok=overwrite)
            if live_view:
                self.enable_live_view()
            if compute_max_intensity_proj:
                self.enable_live_max_intensity_proj_calc()
            self.capture_tiled_image_stack(self.cfg.volume_x_um,
                                           self.cfg.volume_y_um,
                                           self.cfg.volume_z_um,
                                           self.cfg.imaging_wavelengths,
                                           self.cfg.tile_overlap_x_percent,
                                           self.cfg.tile_overlap_y_percent,
                                           self.cfg.tile_prefix,
                                           self.cfg.local_storage_dir,
                                           img_folder)
        finally:
            if live_view:
                self.disable_live_view()
            if compute_max_intensity_proj:
                self.disable_live_max_intensity_proj_calc()
            imaging_log_handler.close()
            self.log.removeHandler(imaging_log_handler)
            # Bundle the log and config files with the dataset.
            self.cfg.save(output_folder, overwrite=overwrite)
            # shutil can't overwrite, so we must delete any prior imaging log
            # in the destination folder if overwriting.
            imaging_log_dest = output_folder/Path(imaging_log_filename)
            if overwrite and imaging_log_dest.exists():
                imaging_log_dest.unlink()
            # We must use shutil because we may be moving files across disks.
            shutil.move(str(imaging_log_filepath), str(output_folder))

    def capture_tiled_image_stack(self, volume_x_um: float, volume_y_um: float,
                                  volume_z_um: float,
                                  wavelengths: list,
                                  tile_overlap_x_percent: float = 15.0,
                                  tile_overlap_y_percent: float = 15.0,
                                  tile_prefix: str = "tile",
                                  local_storage_dir: Path = Path("."),
                                  ext_storage_dir: Path = None):
        """ Capture an array of tiles starting from the current position.
        Note: returns to the starting position even if errors occur
              (except stage-related errors).
        """
        # Iterate through the volume through x, then y, then z.
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
            # Force max intensity worker to timeout if it is running.
            #if max_intensity_worker is not None:
            #   TODO: set some flag here.
            #    max_intensity_worker.join()
            # Normal cleanup.
            self.log.info("Stopping camera(s).")
            for n in self.cfg.ncams:
                self.cam[n].cap_stop()
                self.cam[n].buf_release()
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

                # BUFFERING AND IMAGE CAPTURE FOR DISPIM

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
    #             tif.write(self.cam.buf_getlastframedata(), contiguous=True)
    #             self.image_in_hw_buffer.clear()
    #             images_captured += 1
    #             self.log.debug("Acquired zstack image: "
    #                            f"{images_captured}/{image_count}. "
    #                            f"Read took: {perf_counter() - capture_start:.3f} [s].")

    def setup_imaging_for_laser(self, wavelength: int):
        """Configure system to image with the desired laser wavelength.

        Note: this fn changes the filter wheel too.
        """
        # Bail early if this laser is already setup.
        if self.active_laser == wavelength:
            return
        self.log.info(f"Configuring {wavelength} laser.")
        if self.active_laser is not None:
            self.lasers[self.active_laser].disable()
        fw_index = self.cfg.laser_specs[str(wavelength)]['filter_index']
        self.filter_wheel.set_index(fw_index)
        camera_focus_pos = self.cfg.laser_specs[str(wavelength)]['camera_axis_position']
        self.camera_pose.move_absolute(camera_focus_pos, wait=True)
        # Reprovision the DAQ.
        self.configure_daq(wavelength)
        self.active_laser = wavelength
        self.lasers[self.active_laser].enable()

    def enable_live_view(self):
        """Enable live image viewing in a separate thread."""
        # Bail early if it's already enabled.
        if self.live_view_worker is not None:
            return
        self.log.debug("Enabling live view.")
        self.live_view_worker = Thread(target=self._live_view_worker,
                                       daemon=True)

    def disable_live_view(self, wait: bool = False):
        if self.live_view_worker is not None:
            wait_cond = "" if wait else "not "
            self.log.debug(f"Disabling live view and {wait_cond}waiting.")
            self.live_view_enabled.clear()
        if wait:
            self.live_view_worker.join()

    def _live_view_worker(self):
        """Threadworker for viewing images."""
        raise NotImplementedError

    def enable_live_max_intensity_proj_calc(self):
        raise NotImplementedError

    def disable_live_max_intensity_proj_calc(self):
        raise NotImplementedError

    def _live_max_proj_comp_worker(self):
        """Threadworker for live computation of zstack max intensity projection."""
        # TODO: this needs to happen in a separate process with queued images
        #       so as not to slow down the main thread.
        # TODO: error if we cannot compute max intensity projection
        # fast enough to do it live. Bail but don't kill the run.
        raise NotImplementedError

    def reload_config(self):
        """Reload the toml file."""
        self.cfg.reload()

    def close(self):
        """Safely close all open hardware connections."""
        self.log.info("Closing daq.")
        self.daq.close()
        self.active_laser = None
        if not self.simulated:  # These devices cannot be stubbed out.
            self.log.info("Powering down lasers.")
            for wavelength, laser in self.lasers.items():
                self.log.info(f"Powering down {wavelength}[nm] laser.")
                laser.disable()
            self.log.info("Closing ETL.")
            self.etl.close(soft_close=True)
            self.log.info("Closing camera.")
            self.cam.dev_close()
            self.log.info("De-initializing Dcam API.")
            # Dcamapi.uninit can be called repeatedly but we must call it at least
            # once at the end.
            if not self.simulated:
                Dcamapi.uninit()
        self.log.info("Ending log.")
        for handler in self.log_handlers:
            handler.close()
