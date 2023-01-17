import logging
logging.getLogger().setLevel(logging.DEBUG)
logging.basicConfig(filename=r"C:\Acquire Test\imaging_logging.log", encoding='utf-8',level=logging.DEBUG)

import calliphlox
from calliphlox import DeviceKind, Trigger, SampleType, TriggerEvent, SignalIOKind, TriggerEdge, Direction
from pathlib import Path
from datetime import date
import shutil
from time import sleep
import sys
from datetime import datetime

runtime = calliphlox.Runtime()

class AcquireTest():

    def __init__(self):

        global runtime
        logFile = open(r"C:\Acquire Test\imaging_logging.log", 'a')

        self.initialize_camera()
        self.setup_camera()

        run_number = 0
        while True:

            logging.info("Starting cameras.")
            print('Starting Cameras')
            runtime.start()

            total_frames = 20000
            frames_collected = 0
            while total_frames > frames_collected:

                if a := runtime.get_available_data(0):
                    packet = a.get_frame_count()
                    f = next(a.frames())
                    im = f.data().squeeze().copy()

                    f = None  # <-- fails to get the last frames if this is held?
                    a = None  # <-- fails to get the last frames if this is held?

                    frames_collected += packet
                    print(f'frames collected: {frames_collected}')
            print('Total frames collected')
            logging.info('Total frames collected')

            runtime.stop()
            print('Camera stop')
            logging.info('Camera stop')

            print(f'Completed run {run_number} at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")} ')
            logging.info(f'Completed run {run_number} at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")} ')
            run_number += 1

    def initialize_camera(self):

        print('Initializing Camera')
        logging.info('Initializing Camera')
        dm = runtime.device_manager()
        self.p = runtime.get_configuration()
        self.cameras = [
            d.name
            for d in dm.devices()
            if (d.kind == DeviceKind.Camera) and ("C15440" in d.name)
        ]

    def setup_camera(self):

        dm = runtime.device_manager()

        self.p.video[0].camera.identifier = dm.select(DeviceKind.Camera, self.cameras[0])
        self.p.video[0].storage.identifier = dm.select(DeviceKind.Storage, "Trash")
        self.p.video[0].camera.settings.binning = 1
        self.p.video[0].camera.settings.shape = (2304, 2304)
        self.p.video[0].frame_average_count = 0  # disables
        self.p.video[0].max_frame_count = 10000000000
        print('Setting Configuration')
        logging.info('Setting Configuration')
        runtime.set_configuration(self.p)


if __name__ == "__main__" :

    Test = AcquireTest()
