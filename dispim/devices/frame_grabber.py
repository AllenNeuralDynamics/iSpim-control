import logging
from mock import Mock

try:
    import calliphlox
    from calliphlox import DeviceKind, Trigger, SampleType, TriggerEvent, SignalIOKind, TriggerEdge
except ImportError:
    print("WARNING: failed to import calliphlox")
from pathlib import Path


class FrameGrabber:

    def __init__(self):

        self.runtime = calliphlox.Runtime()
        self.dm = self.runtime.device_manager()
        self.p = self.runtime.get_configuration()
        self.log = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

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
            self.p.video[stream_id].storage.identifier = self.dm.select(DeviceKind.Storage, "Trash")
            self.p.video[stream_id].camera.settings.binning = 1
            self.p.video[stream_id].camera.settings.shape = (tile_shape[1], tile_shape[0])
            self.p.video[stream_id].camera.settings.pixel_type = SampleType.U16
            self.p.video[stream_id].frame_average_count = 0  # disables

    def setup_stack_capture(self, output_path: Path, frame_count: int):
        """Setup capturing for a stack. Including tiff file storage location

        :param tile_shape: size 2 tuple of (columns, rows) for single tile
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
            self.p.video[stream_id].camera.settings.triggers = Trigger(
                enable='True',
                line=0,
                event='AcquisitionStart',
                kind='Input',
                edge='Rising')
        self.runtime.set_configuration(self.p)

    def setup_live(self):
        """Setup for live view. Images are sent to trash and there is no max frame count"""

        print(self.p.video[0].camera.settings.triggers)
        for stream_id in range(0, 2):
            self.p.video[stream_id].storage.identifier = self.dm.select(DeviceKind.Storage, "Trash")
            self.p.video[stream_id].max_frame_count = 1000000
            live_trigger_0 = Trigger(enable='True',
                                     line = 2,
                                     event='FrameStart',
                                     kind='Input',
                                     edge='Rising')
            # External Trigger is index 1 in triggers list. Setup dummy trigger to skip index 0
            self.p.video[stream_id].camera.settings.triggers = [Trigger(),live_trigger_0]
            self.p.video[stream_id].camera.settings.triggers[1].enable = True
            self.p.video[stream_id].camera.settings.triggers[1].event = TriggerEvent.FrameStart
            self.p.video[stream_id].camera.settings.triggers[1].kind = SignalIOKind.Input
            self.p.video[stream_id].camera.settings.triggers[1].edge = TriggerEdge.Rising
        out = self.runtime.set_configuration(self.p)
        from pprint import pprint
        pprint(out.dict())


    def get_exposure_time(self):
        exposure_time = [
            self.p.video[stream_id].camera.settings.exposure_time_us
            for stream_id in range(0, 2)
        ]
        return exposure_time

    def set_exposure_time(self, exp_time: float, live: bool = False):
        for stream_id in range(0, 2):
            self.p.video[stream_id].camera.settings.exposure_time_us = exp_time
            print(f'exposure {stream_id} set to: {self.p.video[stream_id].camera.settings.exposure_time_us}')
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
        for stream_id in range(0, 2):
            self.p.video[stream_id].camera.settings.line_interval_us = line_int
            print(f'line interval {stream_id} set to: {self.p.video[stream_id].camera.settings.line_interval_us}')
        if live:
            self.stop()
            self.runtime.set_configuration(self.p)
            self.start()
        else:
            self.runtime.set_configuration(self.p)
        print(self.runtime.get_configuration())

    def start(self):
        """start the setup frame acquisition."""
        self.log.info("Starting camera.")
        self.runtime.start()

    def stop(self):
        """Stop frame acquisition and file writing."""
        self.log.info("Stopping camera.")
        self.runtime.stop()

    def close(self):
        self.runtime = None
