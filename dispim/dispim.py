"""Abstraction of the DISPIM Instrument."""

import numpy as np
from pathlib import Path
from time import perf_counter
from mock import NonCallableMock as Mock
from dispim.dispim_config import DispimConfig
from dispim.devices.camera import Camera
from tigerasi.tiger_controller import TigerController, UM_TO_STEPS
from tigerasi.sim_tiger_controller import TigerController as SimTiger
# TODO: consolodate these later.
from mesospim.spim_base import Spim
from mesospim.devices.tiger_components import SamplePose


class Dispim(Spim):

    def __init__(self, config_filepath: str, log_filename: str = 'debug.log',
                 console_output: bool = True,
                 color_console_output: bool = False,
                 console_output_level: str = 'info', simulated: bool = False):

        super().__init__(config_filepath, log_filename, console_output,
                         color_console_output, console_output_level, simulated)
        self.cfg = DispimConfig(config_filepath)

        # Separate Processes.

        # Hardware
        #self.cam = Camera() if not self.simulated else Mock(Camera())
        #self.ni = WaveformGenerator()
        self.tigerbox = TigerController(**self.cfg.tiger_obj_kwds) if not \
            self.simulated else SimTiger(**self.cfg.tiger_obj_kwds)
        self.sample_pose = SamplePose(self.tigerbox)

        # Extra Internal State attributes for the current image capture
        # sequence. These really only need to persist for logging purposes.
        self.image_index = 0  # current image to capture.
        self.total_tiles = 0  # tiles to be captured.
        self.stage_x_pos = None
        self.stage_y_pos = None

        # Setup hardware according to the config.
        self._setup_motion_stage()
        self._setup_camera()

    def _setup_motion_stage(self):
        """Configure the sample stage for the Exaspim according to the config."""
        # The tigerbox x axis is the sample pose z axis.
        #   TODO: map this in the config.
        self.tigerbox.set_axis_backlash(x=0)
        # TODO: handle axis remapping here. Remapping should probably come from
        #   the config.
        # Tiger X is Tiling Z, Tiger Y is Tiling X, Tiger Z is Tiling Y
        #self.sample_pose.apply_axis_remapping(self.cfg.axis_map)

    def _setup_camera(self):
        """Configure the camera for the Exaspim according to the config."""
        # pass config parameters into object here.
        # TODO: add calliphlox setup stuff here.
        pass

    def run_from_config(self):
        self.collect_volumetric_image(self.cfg.volume_x_um,
                                      self.cfg.volume_y_um,
                                      self.cfg.volume_z_um,
                                      self.cfg.laser_wavelengths,  # TODO: fix this.
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
        pass

    def livestream(self):
        pass

    def close(self):
        """Safely close all open hardware connections."""
        # stuff here.
        super().close()
