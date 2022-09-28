"""Abstraction of the DISPIM Instrument."""

import logging
import numpy as np
from pathlib import Path
from time import perf_counter
from mock import NonCallableMock as Mock
from dispim.dispim_config import DispimConfig
from dispim.devices.frame_grabber import FrameGrabber
from dispim.devices.ni import WaveformHardware
from dispim.compute_waveforms import generate_waveforms
from tigerasi.tiger_controller import TigerController, UM_TO_STEPS
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

        # Hardware
        self.frame_grabber = FrameGrabber() if not self.simulated else \
            Mock(FrameGrabber)
        self.ni = WaveformHardware() if not self.simulated else \
            Mock(WaveformHardware)
        self.tigerbox = TigerController(**self.cfg.tiger_obj_kwds) if not \
            self.simulated else SimTiger(**self.cfg.tiger_obj_kwds)
        self.sample_pose = SamplePose(self.tigerbox,
                                      **self.cfg.sample_pose_kwds)
        # TODO, setup oxxius laser

        # Extra Internal State attributes for the current image capture
        # sequence. These really only need to persist for logging purposes.
        self.image_index = 0  # current image to capture.
        self.total_tiles = 0  # tiles to be captured.
        self.stage_x_pos = None
        self.stage_y_pos = None

        # Setup hardware according to the config.
        self._setup_camera()
        self._setup_lasers()
        self._setup_motion_stage()
        #self._setup_waveform_hardware()

        # TODO, setup cameras with CPX -> frame_grabber()
    def _setup_camera(self):
            pass

        # TODO, set laser powers and modulation mode for oxxius lasers
    def _setup_lasers(self):
        pass

    def _setup_motion_stage(self):
        """Configure the sample stage for the dispim according to the config."""
        # The tigerbox x axis is the sample pose z axis.
        self.sample_pose.set_axis_backlash(z=0.0)
        # Note: Tiger X is Tiling Z, Tiger Y is Tiling X, Tiger Z is Tiling Y.
        #   This axis remapping is handled upon SamplePose __init__.
        # loop over axes and verify in external mode
        for _, axis in self.cfg.tiger_specs['axes'].items():
            pass
            #self.tigerbox.pm(axis, 1)
        # TODO, think about where to store this mapping in config
        # TODO, merge dispim commands in tigerasi
        # TODO, how to call this? via tigerbox?
        # set card 31 (XY stage), 'X" (input), TTL to value of 1
        #self.tigerbox.ttl(31, 'X', 1)
        # TODO, this needs to be buried somewhere else
        # TODO, how to store card # mappings, in config?

    def _setup_waveform_hardware(self, active_wavelength: int):
        # Compute voltages_t.
        # Write it to hardware.
        self.log.error("Writing waveforms to hardware not yet implemented.")
        self.ni.configure(self.cfg.get_daq_cycle_time(), self.cfg.daq_ao_names_to_channels)
        _, voltages_t = generate_waveforms(self.cfg, active_wavelength)
        self.ni.assign_waveforms(voltages_t)

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

        x_grid_step_um = \
            (1 - tile_overlap_x_percent / 100.0) * self.cfg.tile_size_x_um
        y_grid_step_um = \
            (1 - tile_overlap_y_percent / 100.0) * self.cfg.tile_size_y_um

        # Compute step count.
        # Always round up so that we cover the desired imaging region.
        xsteps = ceil((volume_x_um - self.cfg.tile_size_x_um)
                      / x_grid_step_um)
        ysteps = ceil((volume_y_um - self.cfg.tile_size_y_um)
                      / y_grid_step_um)
        zsteps = ceil((volume_z_um - self.cfg.z_step_size_um)
                      / self.cfg.z_step_size_um)

        xtiles, ytiles, ztiles = (1 + xsteps, 1 + ysteps, 1 + zsteps)
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

        # Set the sample starting location as the origin.
        self.sample_pose.home_in_place()
        # Disable backlash compensation on Z to minimize z axis settling time.
        # Disabling Z compensation is ok since we scan in one direction.
        self.sample_pose.set_axis_backlash(z=0.0)
        # Return to start position and kill the frame grabber.
        # Apply a lead-in-move to take out backlash.
        z_backup_pos = -UM_TO_STEPS * self.cfg.stage_backlash_reset_dist_um
        self.log.debug("Applying extra move to take out backlash.")
        self.sample_pose.move_absolute(z=round(z_backup_pos), wait=True)
        self.sample_pose.move_absolute(z=0, wait=True)
        transfer_process = None  # Reference to external tiff transfer process.
        self.stage_x_pos, self.stage_y_pos = (0, 0)
        self.image_index = 1
        # Setup double buffer to capture images while reading out prior image.
        # TODO: setup nidaq here if not already setup.
        try:
            for j in range(ytiles):
                self.stage_x_pos = 0
                self.sample_pose.move_absolute(y=round(self.stage_y_pos),
                                               wait=True)
                for i in range(xtiles):
                    self.sample_pose.move_absolute(x=round(self.stage_x_pos),
                                                   wait=True)
                    # Setup capture of next Z stack.
                    # Open filters, enable active laser; disable the rest.
                    filename = Path(f"{tile_prefix}_{i}_{j}.tiff")
                    filepath_src = local_storage_dir/filename
                    # TODO: consider making step size a fn parameter instead of
                    #   collected strictly from the config.
                    self.log.debug(f"Collecting tile stack: {filename}.")
                    self._collect_stacked_tiff(ztiles, self.cfg.z_step_size_um,
                                               filepath_src)
                    # Start transferring tiff file to its destination.
                    # Note: Image transfer is faster than image capture, but
                    #   we still wait for prior process to finish.
                    if transfer_process is not None:
                        self.log.info("Waiting for tiff transfer process "
                                      "to complete.")
                        transfer_process.join()
                    if img_storage_dir is not None:
                        self.log.info("Starting transfer process for "
                                      f"{filename}.")
                        filepath_dest = img_storage_dir/filename
                        # TODO, use xcopy transfer for speed
                        transfer_process = TiffTransfer(filepath_src,
                                                        filepath_dest)
                        transfer_process.start()
                    self.stage_x_pos += x_grid_step_um * UM_TO_STEPS
                self.stage_y_pos += y_grid_step_um * UM_TO_STEPS
        finally:
            # Wait for waveform playback to finish so we leave signals in their
            # ending voltage states.
            self.ni.wait_until_done()  # use default timeout of 1[s].
            self.log.info("Stopping camera.")
            if transfer_process is not None:
                self.log.debug("joining zstack transfer process.")
                transfer_process.join()
            self.log.info("Returning to start position.")
            self.sample_pose.move_absolute(x=0, y=0, z=0, wait=True)

    def _collect_stacked_tiff(self, tile_count, tile_spacing_um: float,
                              filepath_src):
        # Setup scan
        self.frame_grabber.setup_stack_capture((self.cfg.column_count_px,
                                                self.cfg.row_count_px),
                                               filepath_src, tile_count)

        self.sample_pose.setup_tile_scan('z', 0, tile_count, tile_spacing_um)
        self.frame_grabber.start()
        try:
            self.sample_pose.start_scan()
            # TODO: this func should block until all the frames are captured.
            nframes = 0
            self.log.error("Skipping stack acquisition from FrameGrabber.")
            #while nframes < 100:
            #    # TODO: this should totally get revisited for efficiency.
            #    if a := self.frame_grabber.runtime.get_available_data():
            #        packet = a.get_frame_count()
            #        f = next(a.frames())
            #        im = f.data().squeeze()
            #        nframes += packet
        finally:
            # Return to start position and kill the frame grabber.
            # Apply a lead-in-move to take out backlash.
            z_backup_pos = -UM_TO_STEPS * self.cfg.stage_backlash_reset_dist_um
            self.log.debug("Applying extra move to take out backlash.")
            self.sample_pose.move_absolute(z=round(z_backup_pos), wait=True)
            self.sample_pose.move_absolute(z=0, wait=True)
            self.frame_grabber.stop()

    def livestream(self):
        pass

    def close(self):
        """Safely close all open hardware connections."""
        # stuff here.
        super().close()
