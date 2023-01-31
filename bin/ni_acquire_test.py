import logging
import calliphlox
from calliphlox import DeviceKind, Trigger, SampleType, TriggerEvent, SignalIOKind, TriggerEdge, Direction
from datetime import datetime
from dispim.dispim_config import DispimConfig
from dispim.devices.ni import WaveformHardware
from dispim.compute_waveforms import generate_waveforms
import time


def initialize_camera(runtime):
    print('Initializing Camera')
    logging.info('Initializing Camera')
    dm = runtime.device_manager()
    p = runtime.get_configuration()

    cameras = [
        d.name
        for d in dm.devices()
        if (d.kind == DeviceKind.Camera) and ("C15440" in d.name)
    ]

    return p, cameras


def setup_camera(runtime, p, cameras):
    dm = runtime.device_manager()
    p.video[0].camera.identifier = dm.select(DeviceKind.Camera, cameras[0])
    p.video[0].camera.settings.exposure_time_us = 0
    p.video[0].storage.settings.filename = fr"C:\Acquire Test_0_Zarr.zarr"
    p.video[0].storage.identifier = dm.select(DeviceKind.Storage, 'Zarr')
    p.video[0].camera.settings.binning = 1
    p.video[0].camera.settings.shape = (2304, 2304)
    p.video[0].frame_average_count = 0  # disables
    # p.video[0].max_frame_count = 75
    live_trigger = Trigger(enable='True',
                           line=2,
                           event='FrameStart',
                           kind='Input',
                           edge='Rising')
    # External Trigger is index 1 in triggers list. Setup dummy trigger to skip index 0
    p.video[0].camera.settings.triggers = [Trigger(), live_trigger]

    print('Setting Configuration')
    logging.info('Setting Configuration')
    runtime.set_configuration(p)
    return runtime, p, dm

# Setup NI card

cfg = DispimConfig(r'C:\Users\Administrator\Projects\dispim-control\examples\config.toml')
ni = WaveformHardware(**cfg.daq_obj_kwds)
ni.configure(cfg.get_daq_cycle_time(), cfg.daq_ao_names_to_channels, live = True)
_, voltages_t = generate_waveforms(cfg, 488)
ni.assign_waveforms(voltages_t)

# Setup Camera
runtime = calliphlox.Runtime()
p, cameras = initialize_camera(runtime)
runtime, p, dm = setup_camera(runtime, p, cameras)

# Setup log
#logging.getLogger().setLevel(logging.DEBUG)
logging.basicConfig(filename=r"C:\Acquire Test\example.log", encoding='utf-8', level=logging.DEBUG)

exposure_times = [10]#,15,20,25]
zarr_types = ["Zarr"]#, "ZarrBlosc1ZstdByteShuffle"]

for time in exposure_times:

    logging.info(f"Setting exposure time to {time}")
    print(f"Setting exposure time to {time}")
    p.video[0].camera.settings.exposure_time_us = time
    logging.info(f"Exposure time set to {p.video[0].camera.settings.exposure_time_us}")
    print(f"Exposure time set to {p.video[0].camera.settings.exposure_time_us}")


    for type in zarr_types:
        logging.info(f"Setting file saving to {type}")
        print(f"Setting file saving to {type}")
        p.video[0].storage.settings.filename = fr"C:\Acquire Test_{time}_{type}.zarr"
        p.video[0].storage.identifier = dm.select(DeviceKind.Storage, type)
        runtime.set_configuration(p)

        logging.info("Starting NI card")
        print('Starting NI card')
        ni.start()

        logging.info("Starting cameras.")
        print('Starting Cameras')
        runtime.start()

        total_frames = 1000
        frames_collected = 0
        no_frames = 0
        while total_frames > frames_collected:

            if a := runtime.get_available_data(0):
                packet = a.get_frame_count()
                f = next(a.frames())
                im = f.data().squeeze().copy()

                f = None  # <-- fails to get the last frames if this is held?
                a = None  # <-- fails to get the last frames if this is held?

                frames_collected += packet
                no_frames = 0
                print(f'frames collected: {frames_collected}')

            # else:
            #     no_frames += 1
            #
            # if no_frames > 10000:
            #     print(f'Stopped collecting frames with exposure time: {time} and file storage: {type}. '
            #                  f'Breaking out of acquisition loop')
            #     logging.info(f'Stopped collecting frames with exposure time: {time} and file storage: {type}. '
            #                  f'Breaking out of acquisition loop')
            #     print(no_frames)
            #     break




        print('Total frames collected')
        logging.info('Total frames collected')

        runtime.abort()
        print('Camera stop')
        logging.info('Camera stop')

        print(f'Completed run with exposure time: {time} and file storage: {type} at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")} ')
        logging.info(f'Completed run exposure time: {time} and file storage: {type} at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")} ')

        logging.info("Stopping NI card")
        print('Stopping NI card')
        ni.stop()

print('Finished')

