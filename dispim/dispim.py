"""Abstraction of the DISPIM Instrument."""

import logging
import numpy as np
from pathlib import Path
from time import perf_counter, sleep
from mock import NonCallableMock as Mock
from dispim.dispim_config import DispimConfig
from dispim.devices.frame_grabber import FrameGrabber
from dispim.devices.ni import WaveformHardware
from dispim.compute_waveforms import generate_waveforms
from tigerasi.tiger_controller import TigerController, UM_TO_STEPS
# from tigerasi.device_codes import ControlMode
from tigerasi.sim_tiger_controller import TigerController as SimTiger
# TODO: consolidate these later.
from mesospim.spim_base import Spim
from mesospim.devices.tiger_components import SamplePose
from math import ceil
from mesospim.tiff_transfer import TiffTransfer


class Dispim(Spim):

    def __init__(self, config_filepath: str,
                 simulated: bool = False):
        #self.log = logging.getLogger(__package__)
        self.log = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        # Log setup is handled in the parent class if we pass in a logger.
        super().__init__(config_filepath, simulated = simulated)
        self.cfg = DispimConfig(config_filepath)

        # Instantiate hardware devices
        self.frame_grabber = FrameGrabber() if not self.simulated else \
            Mock(FrameGrabber)
        self.ni = WaveformHardware(**self.cfg.daq_obj_kwds) if not self.simulated else \
            Mock(WaveformHardware)
        self.tigerbox = TigerController(**self.cfg.tiger_obj_kwds) if not \
            self.simulated else SimTiger(**self.cfg.tiger_obj_kwds)
        self.sample_pose = SamplePose(self.tigerbox,
                                      **self.cfg.sample_pose_kwds)
        # TODO, setup oxxius laser

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

    def _setup_camera(self):
        pass

    def _setup_lasers(self):
        pass

    def _setup_motion_stage(self):
        """Configure the sample stage for the dispim according to the config."""
        self.sample_pose.set_axis_backlash(z=0.0)
        # TODO, set speed of sample Z / tiger X axis to ~1
        # Note: Tiger X is Tiling Z, Tiger Y is Tiling X, Tiger Z is Tiling Y.
        #   This axis remapping is handled upon SamplePose __init__.
        # loop over axes and verify in external mode
        # TODO, think about where to store this mapping in config
        # TODO, merge dispim commands in tigerasi
        # TODO, how to call this? via tigerbox?
        # set card 31 (XY stage), 'X" (input), TTL to value of 1
        # TODO, this needs to be buried somewhere else
        # TODO, how to store card # mappings, in config?

    def _setup_waveform_hardware(self, active_wavelength: int):
        self.log.info("Configuring waveforms for hardware.")
        self.ni.configure(self.cfg.get_daq_cycle_time(), self.cfg.daq_ao_names_to_channels)
        self.log.info("Generating waveforms to hardware.")
        _, voltages_t = generate_waveforms(self.cfg, active_wavelength)
        self.log.info("Writing waveforms to hardware.")
        self.ni.assign_waveforms(voltages_t)
        # TODO: Put all corresponding tigerbox components in external control mode.
        # externally_controlled_axes = \
        #     {a: ControlMode.EXTERNAL_CLOSED_LOOP for a in
        #      self.cfg.ni_controlled_tiger_axes}
        # self.tigerbox.pm(**externally_controlled_axes)

    # TODO: this should be a base class thing.
    def check_ext_disk_space(self, dataset_size):
        self.log.warning("Checking disk space not implemented.")

    def run_from_config(self):
        self.collect_volumetric_image(self.cfg.volume_x_um,
                                      self.cfg.volume_y_um,
                                      self.cfg.volume_z_um,
                                      self.cfg.channels,
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

        # Calculate number of tiles in XYUZ
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

        try:
            for j in range(ytiles):
                # Move X position back to 0
                # Move to specified Y position
                # TODO, set speed of sample Y / tiger Z axis to ~1 mm/s
                self.stage_x_pos = 0
                self.log.info(f"Moving to y={self.stage_y_pos}.")
                self.tigerbox.move_axes_absolute(z=round(self.stage_y_pos), wait_for_output=True, wait_for_reply=True)
                while self.tigerbox.is_moving():
                    pos = self.tigerbox.get_position('Z')
                    self.log.warning(f"Stage is moving! ! Y = {pos['Z']} -> {round(self.stage_y_pos)}")
                    sleep(0.01)

                for i in range(xtiles):
                    # Move to specified X position
                    # TODO, set speed of sample X / tiger Y axis to ~1 mm/s
                    self.log.info(f"Moving to x={self.stage_x_pos}.")
                    self.tigerbox.move_axes_absolute(y=round(self.stage_x_pos), wait_for_output=True, wait_for_reply=True)
                    while self.tigerbox.is_moving():
                        pos = self.tigerbox.get_position('Y')
                        self.log.warning(f"Stage is still moving! X = {pos['Y']} -> {round(self.stage_x_pos)}")
                        sleep(0.01)

                    for ch in channels:
                        # Move to specified Z position
                        # TODO, set speed of sample Z / tiger X axis to ~1 mm/s
                        self.log.info("Applying extra move to take out backlash.")
                        z_backup_pos = -UM_TO_STEPS * self.cfg.stage_backlash_reset_dist_um
                        self.tigerbox.move_axes_absolute(x=round(z_backup_pos), wait_for_output=True, wait_for_reply=True)
                        self.log.info(f"Moving to z={0}.")
                        self.tigerbox.move_axes_absolute(x=0, wait_for_output=True, wait_for_reply=True)
                        while self.tigerbox.is_moving():
                            pos = self.tigerbox.get_position('X')
                            self.log.warning(f"Stage is moving! Z =  {pos['X']} -> {0}")
                            sleep(0.01)
                        # TODO, set speed of sample Z / tiger X axis to ~0.01 mm/s

                        # TODO: setup other channel specific items (filters, lasers)
                        self.log.info(f"Setting up NIDAQ for active channel: {ch}")
                        self._setup_waveform_hardware(ch)

                        # Setup capture of next Z stack.
                        filename = Path(f"{tile_prefix}_X_{i:0>4d}_Y_{j:0>4d}_Z_{0:0>4d}_CH_{ch:0>4d}.tiff")
                        filepath_src = local_storage_dir/filename
                        self.log.info(f"Collecting tile stack: {filepath_src}")

                        # TODO: consider making step size a fn parameter instead of
                        #   collected strictly from the config.
                        tile_position = self.stage_x_pos/UM_TO_STEPS/1000.0 # convert to [mm] units
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
            # TODO, is this needed?
            # self.frame_grabber.close()
            if transfer_process is not None:
                self.log.info("Joining file transfer process.")
                transfer_process.join()

    def _collect_stacked_tiff(self, tile_position, tile_count, tile_spacing_um: float,
                              filepath_src):
        # TODO: set up slow scan axis to take into account sample pose X location
        try:
            self.log.info(f"Configuring framegrabber")
            self.frame_grabber.setup_stack_capture((self.cfg.sensor_column_count,
                                                    self.cfg.sensor_row_count),
                                                   filepath_src, tile_count)
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

    def livestream(self):
        pass

    def close(self):
        """Safely close all open hardware connections."""
        # stuff here.
        super().close()
