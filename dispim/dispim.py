"""Abstraction of the DISPIM Instrument."""

# FIXME: this should live in the napari gui side, and the function
#   we want to launch in a napari thread should exist as a standalone
#   option here.
from napari.qt.threading import thread_worker
import logging
import numpy as np
from pathlib import Path
from time import perf_counter, sleep
from mock import NonCallableMock as Mock
from threading import Thread, Event
from collections import deque
from dispim.dispim_config import DispimConfig
from dispim.devices.frame_grabber import FrameGrabber
from dispim.devices.ni import WaveformHardware
from dispim.compute_waveforms import generate_waveforms
from dispim.devices.oxxius_components import LaserHub
from serial import Serial, EIGHTBITS, STOPBITS_ONE, PARITY_NONE
from tigerasi.tiger_controller import TigerController, UM_TO_STEPS
from tigerasi.sim_tiger_controller import TigerController as SimTiger
# TODO: consolidate these later.
from mesospim.spim_base import Spim
from mesospim.devices.tiger_components import SamplePose
from math import ceil
from mesospim.tiff_transfer import TiffTransfer

OXXIUS_COM_SETUP = \
            {
                "baudrate": 9600,
                "bytesize": EIGHTBITS,
                "parity": PARITY_NONE,
                "stopbits": STOPBITS_ONE,
                "xonxoff": False,
                "timeout": 1
            }

class Dispim(Spim):

    def __init__(self, config_filepath: str,
                 simulated: bool = False):
        self.log = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        # Log setup is handled in the parent class if we pass in a logger.
        super().__init__(config_filepath, simulated=simulated)
        self.cfg = DispimConfig(config_filepath)
        # Instantiate hardware devices
        self.frame_grabber = FrameGrabber() if not self.simulated else \
            Mock(FrameGrabber)
        self.ni = WaveformHardware(**self.cfg.daq_obj_kwds) if not self.simulated else \
            Mock(WaveformHardware)
        self.tigerbox = TigerController(**self.cfg.tiger_obj_kwds) if not \
            self.simulated else SimTiger(**self.cfg.tiger_obj_kwds)
        self.sample_pose = SamplePose(self.tigerbox,)
                                      #**self.cfg.sample_pose_kwds)
        # TODO, setup oxxius laser
        self.lasers = {}  # populated in _setup_lasers.

        # Extra Internal State attributes for the current image capture
        # sequence. These really only need to persist for logging purposes.
        self.total_tiles = 0  # tiles to be captured.
        self.stage_x_pos = None
        self.stage_y_pos = None

        # Setup hardware according to the config.
        self._setup_camera()
        self._setup_lasers()
        self._setup_motion_stage()
        # TODO, setup cameras with CPX -> frame_grabber()
        # TODO, note NIDAQ is channel specific and gets instantiated within imaging loop

        # Internal state attributes.
        self.active_laser = None  # Bookkeep which laser is configured.
        self.live_status = None  # Bookkeep if we are running in live mode.

        self.livestream_worker = None  # captures images during livestream
        self.livestream_enabled = Event()
        self.image_in_hw_buffer = Event()
        self.volumetric_imaging = Event()
        self.img_deque = deque(maxlen=2)  # circular buffer

        # derived constants.
        self.IMG_MIN = 0
        self.IMG_MAX = (1 << (8 * self.cfg.image_dtype.itemsize)) - 1

        # camera id number
        self.stream_id = 0
        self.not_stream_id = 1

    def _setup_camera(self):
        """Configure general settings"""
        self.frame_grabber.setup_cameras((self.cfg.sensor_column_count,
                                          self.cfg.sensor_row_count))

    def _setup_lasers(self):
        """Setup lasers that will be used for imaging. Warm them up, etc."""

        self.ser = Serial(port = 'COM7', **OXXIUS_COM_SETUP) if not self.simulated else None

        for wl, specs in self.cfg.laser_specs.items():
            self.lasers[int(wl)] = LaserHub(specs['prefix'], self.ser) if not self.simulated \
                else Mock(LaserHub)

            self.log.debug(f"Setting up {specs['color']} laser.")

        print(self.lasers)

    def _setup_motion_stage(self):
        """Configure the sample stage for the dispim according to the config."""
        self.log.info("Setting backlash in Z to 0")
        #self.sample_pose.set_axis_backlash(Z=0.0)
        self.log.info("Setting speeds to 1.0 mm/sec")
        #self.tigerbox.set_speed(X=1.0, Y=1.0, Z=1.0)
        # Note: Tiger X is Tiling Z, Tiger Y is Tiling X, Tiger Z is Tiling Y.
        #   This axis remapping is handled upon SamplePose __init__.
        # loop over axes and verify in external mode
        # TODO, think about where to store this mapping in config
        # TODO, merge dispim commands in tigerasi
        # TODO, how to call this? via tigerbox?
        # set card 31 (XY stage), 'X" (input), TTL to value of 1
        # TODO, this needs to be buried somewhere else
        # TODO, how to store card # mappings, in config?

    def _setup_waveform_hardware(self, active_wavelength: int, live: bool = False):

        self.log.info("Configuring waveforms for hardware.")
        self.ni.configure(self.cfg.get_daq_cycle_time(), self.cfg.daq_ao_names_to_channels, live)
        self.log.info("Generating waveforms to hardware.")
        _, voltages_t = generate_waveforms(self.cfg, active_wavelength)
        self.log.info("Writing waveforms to hardware.")
        self.ni.assign_waveforms(voltages_t)
        # TODO: Put all corresponding tigerbox components in external control mode.
        # externally_controlled_axes = \
        #     {a: ControlMode.EXTERNAL_CLOSED_LOOP for a in
        #      self.cfg.ni_controlled_tiger_axes}
        # self.tigerbox.pm(**externally_controlled_axes)

        if self.active_laser is not None and self.live_status == live:
            self.ni.start()

    # TODO: this should be a base class thing.s
    def check_ext_disk_space(self, dataset_size):
        self.log.warning("Checking disk space not implemented.")

    def run_from_config(self):

        if self.livestream_enabled.is_set():
            self.stop_livestream()
        self.collect_volumetric_image(self.cfg.volume_x_um,
                                      self.cfg.volume_y_um,
                                      self.cfg.volume_z_um,
                                      self.cfg.imaging_specs['laser_wavelengths'],
                                      self.cfg.tile_overlap_x_percent,
                                      self.cfg.tile_overlap_y_percent,
                                      self.cfg.tile_prefix,
                                      self.cfg.local_storage_dir,
                                      self.img_storage_dir,
                                      self.deriv_storage_dir)

    def collect_volumetric_image(self, volume_x_um: float, volume_y_um: float,
                                 volume_z_um: float,
                                 channels: list,
                                 tile_overlap_x_percent: float,
                                 tile_overlap_y_percent: float,
                                 tile_prefix: str,
                                 local_storage_dir: Path = Path("."),
                                 img_storage_dir: Path = None,
                                 deriv_storage_dir: Path = None):
        """Collect a tiled volumetric image with specified size/overlap specs.
        """
        # Calculate tiling offsets in XY (we scan along Z)
        x_grid_step_um = \
            (1 - tile_overlap_x_percent / 100.0) * self.cfg.tile_size_x_um
        y_grid_step_um = \
            (1 - tile_overlap_y_percent / 100.0) * self.cfg.tile_size_y_um

        # Calculate number of tiles in XYZ
        # Always round up so that we cover the desired imaging region.
        xtiles = ceil((volume_x_um - self.cfg.tile_size_x_um)
                      / x_grid_step_um)
        ytiles = ceil((volume_y_um - self.cfg.tile_size_y_um)
                      / y_grid_step_um)
        ztiles = ceil((volume_z_um - self.cfg.z_step_size_um)
                      / self.cfg.z_step_size_um)
        xtiles, ytiles, ztiles = (1 + xtiles, 1 + ytiles, 1 + ztiles)
        self.total_tiles = xtiles * ytiles * ztiles * len(channels)

        # Check if we will violate storage directory limits with our
        #   run. (Check local and external storage.)
        gigabytes_per_image = self.cfg.bytes_per_image / (1024 ** 3)
        dataset_gigabytes = self.total_tiles * gigabytes_per_image
        try:  # Log errors if they occur.
            self.check_ext_disk_space(dataset_gigabytes)
        except AssertionError as e:
            self.log.error(e)
            raise
        # TODO, test network speeds?
        # TODO, check that networked storage is visible?

        # TODO, check if stage homing is necessary?
        # Set the sample starting location as the origin.
        # self.sample_pose.home_in_place()

        transfer_process = None  # Reference to external tiff transfer process.
        self.stage_x_pos, self.stage_y_pos = (0, 0)
        self.volumetric_imaging.set()  # TODO: do we need this?
        try:
            for j in range(ytiles):
                # Move X position back to 0
                # Move to specified Y position
                self.log.info("Setting speed in Y to 1.0 mm/sec")
                # TODO: handle this through sample pose class, which remaps axes
                self.tigerbox.set_speed(Z=1.0)
                self.stage_x_pos = 0
                self.log.info(f"Moving to y={self.stage_y_pos}.")
                self.tigerbox.move_axes_absolute(z=round(self.stage_y_pos), wait_for_output=True, wait_for_reply=True)
                while self.tigerbox.is_moving():
                    pos = self.tigerbox.get_position('Z')
                    distance = abs(pos['Z'] - round(self.stage_y_pos))
                    if distance < 1.0:
                        self.tigerbox.halt()
                        break
                    else:
                        self.log.warning(f"Stage is moving! ! Y = {pos['Z']} -> {round(self.stage_y_pos)}")
                        sleep(0.1)

                for i in range(xtiles):
                    # Move to specified X position
                    self.log.info("Setting speed in X to 1.0 mm/sec")
                    # TODO: handle this through sample pose class, which remaps axes
                    self.tigerbox.set_speed(Y=1.0)
                    self.log.info(f"Moving to x={round(self.stage_x_pos)}.")
                    self.tigerbox.move_axes_absolute(y=round(self.stage_x_pos),
                                                     wait_for_output=True,
                                                     wait_for_reply=True)
                    while self.tigerbox.is_moving():
                        pos = self.tigerbox.get_position('Y')
                        distance = abs(pos['Y'] - round(self.stage_x_pos))
                        if distance < 1.0:
                            self.tigerbox.halt()
                            break
                        else:
                            self.log.warning(f"Stage is still moving! X = {pos['Y']} -> {round(self.stage_x_pos)}")
                            sleep(0.1)

                    for ch in channels:

                        self.setup_imaging_for_laser(ch)
                        # Move to specified Z position
                        self.log.info("Setting speed in Z to 1.0 mm/sec")
                        # TODO: handle this through sample pose class, which remaps axes
                        self.tigerbox.set_speed(X=1.0)
                        self.log.info("Applying extra move to take out backlash.")
                        z_backup_pos = -UM_TO_STEPS * self.cfg.stage_backlash_reset_dist_um
                        self.tigerbox.move_axes_absolute(x=round(z_backup_pos),
                                                         wait_for_output=True,
                                                         wait_for_reply=True)
                        self.log.info(f"Moving to z={0}.")
                        self.tigerbox.move_axes_absolute(x=0, wait_for_output=True, wait_for_reply=True)
                        while self.tigerbox.is_moving():
                            pos = self.tigerbox.get_position('X')
                            distance = abs(pos['X'] - round(0))
                            if distance < 1.0:
                                self.tigerbox.halt()
                                break
                            else:
                                self.log.warning(f"Stage is moving! Z =  {pos['X']} -> {0}")
                                sleep(0.1)

                        self.log.info("Setting speed in Z to 0.01 mm/sec")
                        self.tigerbox.set_speed(X=0.01)

                        # TODO: setup other channel specific items (filters, lasers)
                        self.log.info(f"Setting up NIDAQ for active channel: {ch}")
                        self._setup_waveform_hardware(ch)

                        # Setup capture of next Z stack.
                        filename = Path(f"{tile_prefix}_X_{i:0>4d}_Y_{j:0>4d}_Z_{0:0>4d}_CH_{ch:0>4d}.tiff")
                        filepath_src = local_storage_dir/filename
                        self.log.info(f"Collecting tile stack: {filepath_src}")

                        # TODO: consider making step size a fn parameter instead of
                        #   collected strictly from the config.
                        tile_position = self.stage_x_pos/UM_TO_STEPS/1000.0  # convert to [mm] units
                        self.log.info(f"Tile position [mm]: {tile_position}.")
                        self.log.info(f"Tiles [#]: {ztiles}.")
                        self.log.info(f"Tile spacing [um]: {self.cfg.z_step_size_um}.")
                        self._collect_stacked_tiff(tile_position, ztiles, self.cfg.z_step_size_um,
                                                   filepath_src)

                        # Start transferring tiff file to its destination.
                        # Note: Image transfer is faster than image capture, but
                        #   we still wait for prior process to finish.
                        if transfer_process is not None:
                            self.log.info("Waiting for tiff transfer process "
                                          "to complete.")
                            transfer_process.join()
                        if img_storage_dir is not None:
                            filepath_dest = img_storage_dir/filename
                            self.log.info("Starting transfer process for "
                                          f"{filepath_dest}.")
                            # TODO, use xcopy transfer for speed
                            transfer_process = TiffTransfer(filepath_src,
                                                            filepath_dest)
                            transfer_process.start()

                        # TODO, set speed of sample Z / tiger X axis to ~1

                    self.stage_x_pos += x_grid_step_um * UM_TO_STEPS
                self.stage_y_pos += y_grid_step_um * UM_TO_STEPS

        finally:
            # self.log.info("Returning to start position.")
            # self.sample_pose.move_absolute(x=0, y=0, z=0, wait=True)
            self.log.info(f"Closing camera")
            self.frame_grabber.close()
            self.volumetric_imaging.clear()  # TODO: Do we need this?
            if transfer_process is not None:
                self.log.info("Joining file transfer process.")
                transfer_process.join()

    def _collect_stacked_tiff(self, tile_position, tile_count, tile_spacing_um: float,
                              filepath_src):
        # TODO: set up slow scan axis to take into account sample pose X location
        try:
            self.log.info(f"Configuring framegrabber")
            self.frame_grabber.setup_stack_capture(filepath_src, tile_count)
            self.log.info(f"Configuring stage scan parameters")
            self.sample_pose.setup_tile_scan('z', 0, tile_count, tile_spacing_um, tile_position)

            self.log.info(f"Starting framegrabber")
            self.frame_grabber.start()
            self.log.info(f"Starting NIDAQ")
            self.ni.start()
            self.log.info(f"Starting stage")
            self.sample_pose.start_scan()

            frames = 0
            while frames < tile_count:
                frames = self.ni.counter_task.read()
                self.log.info(f"{frames} frames captured")
                sleep(0.01)

            sleep(10)

        finally:
            self.log.info(f"Stopping NIDAQ")
            self.ni.stop()
            self.log.info(f"Closing NIDAQ")
            self.ni.close()
            self.log.info(f"Stopping framegrabber")
            self.frame_grabber.stop()

    def start_livestream(self, wavelength: int):
        """Repeatedly play the daq waveforms and buffer incoming images."""
        if wavelength not in self.cfg.imaging_specs['laser_wavelengths']:
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
        self.frame_grabber.setup_live()
        self.ni.start()
        self.livestream_enabled.set()

    def stop_livestream(self, wait: bool = False):
        # Bail early if it's already stopped.
        if not self.livestream_enabled.is_set():
            self.log.warning("Not starting. Livestream is already stopped.")
            return
        wait_cond = "" if wait else "not "
        self.log.debug(f"Disabling livestream and {wait_cond}waiting.")
        self.livestream_enabled.clear()
        self.frame_grabber.stop()

        self.ni.stop()
        self.ni.close()
        self.live_status = False  # TODO: can we get rid of this if we're always stopping the livestream?
        self.lasers[self.active_laser].disable()
        self.active_laser = None


    @thread_worker
    def _livestream_worker(self):
        """Pulls images from the camera and puts them into the ring buffer."""
        image_wait_time = round(5 * self.cfg.get_daq_cycle_time() * 1e3)
        self.frame_grabber.start()  # ?
        while self.livestream_enabled.is_set():
            self.stream_id, self.not_stream_id = self.not_stream_id, self.stream_id
            if self.simulated:
                sleep(1 / 16)
                blank = np.zeros((self.cfg.sensor_row_count,
                                  self.cfg.sensor_column_count),
                                 dtype=self.cfg.image_dtype)
                noise = np.random.normal(0, .1, blank.shape)
                yield noise + blank, self.stream_id

            elif packet := self.frame_grabber.runtime.get_available_data(self.stream_id):
                f = next(packet.frames())
                im = f.data().squeeze().copy()  # TODO: copy?

                f = None  # <-- will fail to get the last frames if this is held?
                packet = None  # <-- will fail to get the last frames if this is held?

                yield im, self.stream_id


    def setup_imaging_for_laser(self, wavelength: int, live: bool = False):
        """Configure system to image with the desired laser wavelength.
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
        # Reprovision the DAQ.
        self._setup_waveform_hardware(wavelength, live)
        self.active_laser = wavelength
        self.lasers[self.active_laser].enable()

    def move_sample_absolute(self, x: int = None, y: int = None, z: int = None):
         """Convenience function for moving the sample from a UI."""
         self.sample_pose.move_absolute(x=x, y=y, z=z, wait=True)

    def move_sample_relative(self, x: int = None, y: int = None, z: int = None):
        """Convenience func for moving the sample from a UI (units: steps)."""
        self.sample_pose.move_relative(x=x, y=y, z=z, wait=True)

    def get_sample_position(self):
        return self.sample_pose.get_position()

    def close(self):
        """Safely close all open hardware connections."""
        # stuff here.
        self.ni.close()
        for wavelength, laser in self.lasers.items():
            self.log.info(f"Powering down {wavelength}[nm] laser.")
            laser.disable()
        self.ser.close()  # TODO: refactor oxxius lasers into wrapper class.
        self.tigerbox.ser.close()
        self.frame_grabber.close()
        super().close()
