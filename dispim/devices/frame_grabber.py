import logging
from mock import Mock

try:
    import calliphlox
    from calliphlox import DeviceKind, Trigger, SampleType
except ImportError:
    print("WARNING: failed to import calliphlox")
from pathlib import Path


class FrameGrabber:

    def __init__(self):
        self.runtime = calliphlox.Runtime()
        self.dm = self.runtime.device_manager()
        self.p = self.runtime.get_configuration()

        self.cameras = [
            d.name
            for d in self.dm.devices()
            if (d.kind == DeviceKind.Camera) and ("C15440" in d.name)
        ]

        self.log = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

    def setup_cameras(self, tile_shape: tuple):
        """General setup for both cameras for both livestream and stack capture

        :param tile_shape: size 2 tuple of (columns, rows) for single tile"""

        for stream_id in range(0, 2):
            self.p.video[stream_id].camera.identifier = self.dm.select(DeviceKind.Camera, self.cameras[stream_id])
            self.p.video[stream_id].camera.settings.binning = 1
            self.p.video[stream_id].camera.settings.shape = (tile_shape[1], tile_shape[0])
            self.p.video[stream_id].camera.settings.pixel_type = SampleType.U16
            self.p.video[stream_id].frame_average_count = 0  # disables

    def setup_stack_capture(self, output_path: Path, frame_count: int):
        """Setup capturing for a stack. Including tiff file storage location

        :param output_path: path where tiff will be saved
        :param frame_count: how many tiles to grab from camera
        """
        # TODO: Should this be looped over so we can configure both cameras at the same time?
        # is there ever a time where there would be different configurations for stack capture?

        for stream_id in range(0, 2):
            self.log.info(f"Configuring camera {stream_id}.")
            self.p.video[stream_id].storage.identifier = self.dm.select(DeviceKind.Storage, "Tiff")
            self.log.info(str(output_path.absolute()))
            self.p.video[stream_id].storage.settings.filename = str(output_path.absolute())
            self.p.video[stream_id].max_frame_count = frame_count

        self.runtime.set_configuration(self.p)

    def setup_live(self):
        """Setup for live view. Images are sent to trash and there is no max frame count"""

        for stream_id in range(0, 2):
            self.p.video[stream_id].storage.identifier = self.dm.select(DeviceKind.Storage, "Trash")
            # self.p.video[stream_id].max_frame_count = inf
        self.runtime.set_configuration(self.p)

    @property
    def exposure_time(self, stream_id: int):
        return self.p.video[stream_id].camera.settings.exposure_time_us

    @exposure_time.setter
    def exposure_time(self, stream_id: int, exp_time: float):
        self.p.video[stream_id].camera.settings.exposure_time_us = exp_time

    @property
    def line_interval(self, stream_id: int):
        return self.p.video[stream_id].camera.settings.line_interval_us

    @line_interval.setter
    def line_interval(self, stream_id: int, line_int: float):
        self.p.video[stream_id].camera.settings.line_interval_us = line_int

    def start(self):
        """start the setup frame acquisition."""
        self.log.info("Starting camera.")
        self.runtime.start()

    def stop(self):
        """Stop frame acquisition and file writing."""
        self.log.info("Stopping camera.")
        self.runtime.stop()

    def close(self):
        pass
