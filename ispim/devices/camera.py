import logging
from acquire import DeviceKind, Trigger, SampleType, Trigger, SignalIOKind, TriggerEdge, Direction, Runtime,\
    AvailableData
from pathlib import Path
from spim_core.devices import camera_core

class Camera(camera_core):

    def __init__(self, name):
        self.log = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self.runtime = Runtime()
        dm = self.runtime.device_manager()
        self.p = self.runtime.get_configuration()
        self.device = None
        for d in dm.devices():
            if (d.kind == DeviceKind.Camera) and (name in d.name):
                self.device = d
                break
        if self.device == None:
            self.log.error(f"Cannot find camera with the name {name}")
            raise

    def configure(self, tile_shape: list, scan_direction: str, exposure_time_ms: float, slit_width_px: int):
        """Configure camera based on values in config"""

        dm = self.runtime.device_manager()
        self.p.video[0].camera.identifier = dm.select(DeviceKind.Camera, self.device)
        self.p.video[0].storage.identifier = dm.select(DeviceKind.Storage, "Trash")
        self.p.video[0].camera.settings.binning = 1
        self.p.video[0].camera.settings.line_interval_us = (exposure_time_ms * 1000000) / tile_shape[0]
        self.p.video[0].camera.settings.shape = (tile_shape[0], tile_shape[1])
        self.runtime.set_configuration(self.p)  # Set shape before setting offset
        self.p.video[0].camera.settings.offset = (int((2304 - tile_shape[0]) / 2), int((2304 - tile_shape[1]) / 2))
        self.p.video[0].camera.settings.pixel_type = SampleType.U16
        self.p.video[0].frame_average_count = 0  # disables
        self.p.video[0].camera.settings.readout_direction = Direction.Forward if scan_direction == 'FORWARD' else Direction.Backward
        self.p.video[0].camera.settings.exposure_time_us = slit_width_px * self.p.video[0].camera.settings.line_interval_us

        self.runtime.set_configuration(self.p)

    def start(self):
        """start the setup frame acquisition."""
        self.log.debug("Starting camera.")
        self.runtime.start()

    def stop(self):
        """Stop frame acquisition and file writing."""
        self.log.debug("Stopping camera.")
        self.runtime.abort()

    def grab_frame(self):
        """Grab latest frame from camera"""
        #TODO check if there's an easier way to get avialable data
        self.frame = AvailableData.frames()
        return self.frame

    def grab_frame_count(self):
        """Grab frame count off camera"""
        return AvailableData.get_frame_count()

    def setup_stack_capture(self, output_path: Path, frame_count: int, filetype: str):
        """Setup capturing for a stack. Including tiff file storage location

        :param output_path: path where camera frame will be saved
        :param frame_count: how many tiles to grab from camera

        """

        dm = self.runtime.device_manager()
        self.log.info(f"Configuring camera for stack capture")
        self.p.video[0].storage.identifier = dm.select(DeviceKind.Storage, filetype) #zarr compression name = ZarrBlosc1ZstdByteShuffle
        self.p.video[0].storage.settings.filename = str(output_path.absolute())
        self.p.video[0].max_frame_count = frame_count
        self.p.video[0].camera.settings.input_triggers.frame_start = Trigger(enable=True, line=0, edge="Rising")

        self.runtime.set_configuration(self.p)

    @property
    def exposure_time_us(self):
        """Get exposure time of the camera"""
        return self.p.video[0].camera.settings.exposure_time_us

    @exposure_time_us.setter
    def exposure_time_us(self, time):
        """Set exposure time of the camera"""
        self.p.video[0].camera.settings.exposure_time_us = time

    @property
    def line_interval_us(self):
        """Get line interval of the camera"""
        return self.p.video[0].camera.settings.line_interval_us

    @line_interval_us.setter
    def line_interval_us(self, time):
        """Set line interval of the camera"""
        self.p.video[0].camera.settings.line_interval_us = time

    @property
    def scan_direction(self):
        """Get scanning direction of camera"""
        return self.p.video[0].camera.settings.readout_direction

    @scan_direction.setter
    def scan_direction(self, direction):
        """Get scanning direction of camera"""
        scan_direction = Direction.Forward if direction == 'FORWARD' else Direction.Backward
        self.p.video[0].camera.settings.readout_direction = scan_direction

    def close(self):
        self.runtime = None

