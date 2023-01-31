import logging
from mock import Mock
from pprint import pprint

try:
    import calliphlox
    from calliphlox import DeviceKind, Trigger, SampleType, TriggerEvent, SignalIOKind, TriggerEdge, Direction
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

        for video, camera in zip(self.p.video, self.cameras):
            video.camera.identifier = self.dm.select(DeviceKind.Camera, camera)
            video.storage.identifier = self.dm.select(DeviceKind.Storage, "Trash")
            video.camera.settings.binning = 1
            video.camera.settings.shape = (tile_shape[1], tile_shape[0])
            video.camera.settings.pixel_type = SampleType.U16
            video.frame_average_count = 0  # disables

        # self.p.video[0].camera.settings.readout_direction = Direction.Forward
        # self.p.video[1].camera.settings.readout_direction = Direction.Backward


    def setup_stack_capture(self, output_paths: list[Path], frame_count: int):
        """Setup capturing for a stack. Including tiff file storage location

        :param output_paths: size 2 list of paths where each camera tiff will be saved
        :param frame_count: how many tiles to grab from camera

        """
        # TODO: Should this be looped over so we can configure both cameras at the same time?
        # is there ever a time where there would be different configurations for stack capture?

        for video, path in zip(self.p.video, output_paths):
            self.log.info(f"Configuring camera.")
            video.storage.identifier = self.dm.select(DeviceKind.Storage, "Zarr") #zarr compression name = ZarrBlosc1ZstdByteShuffle
            self.log.info(str(path.absolute()))
            video.storage.settings.filename = str(path.absolute())
            video.max_frame_count = frame_count
            acq_trigger = Trigger(enable='True',
                                     line=2,
                                     event='FrameStart',
                                     kind='Input',
                                     edge='Rising')
            # External Trigger is index 1 in triggers list. Setup dummy trigger to skip index 0
            video.camera.settings.triggers = [Trigger(), acq_trigger]
        self.runtime.set_configuration(self.p)

    def setup_live(self):
        """Setup for live view. Images are sent to trash and there is no max frame count"""

        stream_id = 0
        for video in self.p.video:
            video.storage.identifier = self.dm.select(DeviceKind.Storage, "Trash")
            video.max_frame_count = 1000000
            live_trigger = Trigger(enable='True',
                                     line = 2,
                                     event='FrameStart',
                                     kind='Input',
                                     edge='Rising')
            # External Trigger is index 1 in triggers list. Setup dummy trigger to skip index 0
            video.camera.settings.triggers = [Trigger(),live_trigger]
            self.runtime.set_configuration(self.p)

    def get_exposure_time(self):
        exposure_time = [
            self.p.video[stream_id].camera.settings.exposure_time_us
            for stream_id in range(0, 2)
        ]
        return exposure_time

    def set_exposure_time(self, exp_time: float, live: bool = False):
        for video in self.p.video:
            video.camera.settings.exposure_time_us = exp_time
            print(f'exposure set to: {video.camera.settings.exposure_time_us}')
        if live:
            self.stop()
            self.runtime.set_configuration(self.p)
            self.start()
        else:
            self.runtime.set_configuration(self.p)
        print(self.runtime.get_configuration())

    def get_line_interval(self):
        line_interval = [
            self.p.video[stream_id].camera.settings.line_interval_us
            for stream_id in range(0, 2)
        ]
        return line_interval

    def set_line_interval(self, line_int: float, live: bool = False):
        for video in self.p.video:
            video.camera.settings.line_interval_us = line_int
            print(f'line interval set to: {video.camera.settings.line_interval_us}')
        if live:
            self.stop()
            self.runtime.set_configuration(self.p)
            self.start()
        else:
            self.runtime.set_configuration(self.p)
        print(self.runtime.get_configuration())

    def get_scan_direction(self, stream_id):
        return self.p.video[stream_id].camera.settings.readout_direction

    def set_scan_direction(self, stream_id, direction:str, live: bool = False):
        direction = Direction.Forward if direction == 'FORWARD' else Direction.Backward
        self.p.video[stream_id].camera.settings.readout_direction = direction
        if live:
            self.stop()
            self.runtime.set_configuration(self.p)
            self.start()
        else:
            self.runtime.set_configuration(self.p)

    def start(self):
        """start the setup frame acquisition."""
        self.log.debug("Starting cameras.")
        self.runtime.start()

    def stop(self):
        """Stop frame acquisition and file writing."""
        self.log.debug("Stopping cameras.")
        self.runtime.stop()

    def close(self):
        self.runtime = None
