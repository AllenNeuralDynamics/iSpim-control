import logging
from mock import Mock

try:
    import calliphlox
    from calliphlox import DeviceKind, Trigger
except ImportError:
    print("WARNING: failed to import calliphlox")
from pathlib import Path


class FrameGrabber:

    def __init__(self):
        self.log = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self.runtime = calliphlox.Runtime()

    def setup_stack_capture(self, tile_shape: tuple, output_path: Path, frame_count: int):
        """Setup capturing for a stack. Including tiff file storage location

        :param tile_shape: size 2 tuple of (columns, rows) for single tile
        :param output_path: path where tiff will be saved
        :param frame_count: how many tiles to grab from camera

        """
        self.log.info("Configuring camera.")
        self.dm = self.runtime.device_manager()
        self.p = self.runtime.get_configuration()
        self.p.camera.identifier = self.dm.select(DeviceKind.Camera, name='C15440-20UP')
        self.p.camera.settings.binning = 1
        self.p.camera.settings.shape = (tile_shape[1], tile_shape[0])
        self.p.storage.identifier = self.dm.select(DeviceKind.Storage, name='Tiff')
        self.p.storage.settings.filename = str(output_path.absolute())
        self.p.max_frame_count = frame_count
        self.p.frame_average_count = 0  # disables
        self.p = self.runtime.set_configuration(self.p)

    def start(self):
        """start the setup frame acquisition."""
        self.log.info("Starting camera.")
        self.runtime.start()

    def stop(self):
        """Stop frame acquisition and file writing."""
        self.log.info("Stopping camera.")
        self.runtime.stop()
