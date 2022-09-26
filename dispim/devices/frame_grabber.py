import logging
from mock import Mock
try:
    from calliphlox import DeviceKind, Trigger, Runtime
except ImportError:
    print("WARNING: failed to import calliphlox")
from pathlib import Path


class FrameGrabber:

    def _init_(self):
        self.runtime = Runtime()  # do we need one per stack or one ever?
        self.dm = self.runtime.device_manager()
        self.log = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

    def setup_stack_capture(self, tile_shape: tuple, output_path: Path, frame_count: int):
        """Setup capturing for a stack. Including tiff file storage location

        :param tile_shape: size 2 tuple of (columns, rows) for single tile
        :param output_path: path where tiff will be saved
        :param frame_count: how many tiles to grab from camera

        """
        self.log.error("Skipping FrameGrabber stack capture setup!")
        # # TODO: can we overwrite the runtime configuration over and over, or
        # #   do we need to create a new Runtime object?
        # p = self.runtime.get_configuration()

        # p.camera.identifier = self.dm.select(DeviceKind.Camera,
        #                                      name='C15440-20UP')
        # p.camera.settings.binning = 1
        # p.camera.settings.shape = tile_shape
        # p.storage.identifier = self.dm.select(DeviceKind.Storage, name='Tiff')
        # p.storage.settings.filename = str(output_path.absolute())
        # p.max_frame_count = frame_count
        # p.frame_average_count = 0  # disables
        # #TODO: what is this?
        # self.runtime.set_configuration(p)

    def start(self):
        """start the setup frame acquisition."""
        self.log.error("Skipping FrameGrabber start.")
        #self.runtime.start()
        # TODO: figure out some sensible way to know when this operation is
        #  done?

    def stop(self):
        """Stop frame acquisition and file writing."""
        # TODO: what happens if we call this before capturing the correct
        #  frame count?
        self.log.error("Skipping FrameGrabber stop.")
        #self.runtime.stop()