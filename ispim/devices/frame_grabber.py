import logging
from mock import Mock
from pprint import pprint
import numpy as np
from pathlib import Path
from time import time
from exa_spim_refactor.devices.camera.hamamatsu_dcam import Camera

class FrameGrabber:

    def __init__(self):

        self.camera = Camera("500749")

        self.log = logging.getLogger(f"{__name__}.{self.__class__.__name__}")

    def setup_cameras(self, tile_shape: tuple):
        """General setup for both cameras for both livestream and stack capture

        param tile_shape: size 2 tuple of (columns, rows) for single tile"""

        self.camera.roi = {"width_px": tile_shape[0], "height_px": tile_shape[1]}
        print(self.camera.roi)
        # dm = self.runtime.device_manager()
        # for stream_id, camera in zip(range(0, len(self.cameras)), self.cameras):
        #     self.p.video[stream_id].camera.identifier = dm.select(DeviceKind.Camera, camera)
        #     self.p.video[stream_id].storage.identifier = dm.select(DeviceKind.Storage, "Trash")
        #     self.p.video[stream_id].camera.settings.binning = 1
        #     self.p.video[stream_id].camera.settings.shape = (tile_shape[0], tile_shape[1])
        #     self.runtime.set_configuration(self.p)
        #     self.p.video[stream_id].camera.settings.offset = (int((2304 - tile_shape[0])/2), int((2304 - tile_shape[1])/2)) #TODO: Not hard code 2304
        #     self.p.video[stream_id].camera.settings.pixel_type = SampleType.U16
        #     self.p.video[stream_id].frame_average_count = 0  # disables
        #     self.runtime.set_configuration(self.p)

    def setup_stack_capture(self, output_paths: list[Path], frame_count: int, filetype: str):
        """Setup capturing for a stack. Including tiff file storage location

        :param output_paths: size 2 list of paths where each camera tiff will be saved
        :param frame_count: how many tiles to grab from camera

        """
        self.camera.prepare()
        # dm = self.runtime.device_manager()
        # for stream_id in range(0, len(self.cameras)):
        #     self.log.info(f"Configuring camera.")
        #     self.p.video[stream_id].storage.identifier = dm.select(DeviceKind.Storage, filetype) #zarr compression name = ZarrBlosc1ZstdByteShuffle
        #     self.p.video[stream_id].storage.settings.filename = str(output_paths[stream_id].absolute())
        #     self.p.video[stream_id].max_frame_count = frame_count
        #     self.p.video[stream_id].camera.settings.input_triggers.frame_start = acquire.Trigger(enable=True, line=0, edge="Rising")
        # self.runtime.set_configuration(self.p)
        # print(self.runtime.get_configuration())


    def collect_background(self, frame_average=1):
        """Retrieve a background image as a 2D numpy array with shape (rows, cols). """

        dm = self.runtime.device_manager()
        filetype_id = self.p.video[0].storage.identifier
        self.p.video[0].camera.settings.input_triggers.frame_start = acquire.Trigger(enable=False, line=0,
                                                                                     edge="NotApplicable")
        self.p.video[0].storage.identifier = dm.select(DeviceKind.Storage, 'Trash')
        self.runtime.set_configuration(self.p)
        # Initialize background image array
        bkg_image = np.zeros((frame_average,
                              self.p.video[0].camera.settings.shape[1],
                              self.p.video[0].camera.settings.shape[0]),
                             dtype='uint16')
        # Rapid fire camera and pull desired frame number out
        self.start()
        i = 0
        start_time = time()
        while i < frame_average and time()-start_time<60:
            with self.runtime.get_available_data(0) as a:
                packet = a.get_frame_count()
                if packet == 0:
                    continue
                f = next(a.frames())
                bkg_image[i] = f.data().squeeze().copy()
                i += 1
        self.log.info(f"Averaging {frame_average} background images")
        self.stop()

        self.p.video[0].storage.identifier = filetype_id
        self.p.video[0].camera.settings.input_triggers.frame_start = acquire.Trigger(enable=True, line=0, edge="Rising")
        self.runtime.set_configuration(self.p)
        # Return median averaged 2D background image
        return np.median(bkg_image, axis=0).astype('uint16')

    def get_exposure_time(self):
        return self.camera.exposure_time_ms

    def set_exposure_time(self, exp_time: float, live: bool = False):
        # for video in self.p.video:
        #     video.camera.settings.exposure_time_us = exp_time
        #     self.log.debug(f'exposure set to: {video.camera.settings.exposure_time_us}')
        #
        # if live:
        #     self.stop()
        #     self.runtime.set_configuration(self.p)
        #     self.start()
        # else:
        #     self.runtime.set_configuration(self.p)
        self.camera.exposure_time_ms = exp_time

    def grab_frame(self):
        return self.camera.grab_frame()

    def get_line_interval(self):
        return self.camera.line_interval_us

    def set_line_interval(self, line_int: float, live: bool = False):
        self.camera.line_interval_us = line_int

    def get_scan_direction(self, stream_id):
        direction = self.camera.readout_mode
        READOUT_MODES = {
            "FORWARD":"light sheet forward",
            "BACKWARD":"light sheet backward"
        }
        return READOUT_MODES[direction]

    def set_scan_direction(self, stream_id, direction:str, live: bool = False):
        direction =  "light sheet forward" if direction == 'FORWARD' else "light sheet backward"
        self.camera.readout_mode = direction

    def start(self):
        """start the setup frame acquisition."""
        self.log.debug("Starting cameras.")
        self.camera.start()

    def stop(self):
        """Stop frame acquisition and file writing."""
        self.log.debug("Stopping cameras.")
        self.camera.stop()

    def close(self):
        self.camera.close()
