import logging
import calliphlox
from calliphlox import DeviceKind, Trigger, SampleType, TriggerEvent, SignalIOKind, TriggerEdge, Direction
from datetime import datetime
from dispim.dispim_config import DispimConfig
from dispim.devices.ni import WaveformHardware
from dispim.compute_waveforms import generate_waveforms
from time import time as now
import os

logging.getLogger().handlers.clear()

def initialize_camera(runtime):

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
    p.video[0].camera.settings.line_interval_us = 25
    p.video[0].storage.settings.filename = fr"C:\Acquire Test_0_Zarr.zarr"  # Need to put something at the beginning but can change later
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


    logging.info('Setting Configuration')
    runtime.set_configuration(p)
    return runtime, p, dm

# Setup NI card

cfg = DispimConfig(r'C:\Users\Administrator\Projects\dispim-control\examples\config.toml')
ni = WaveformHardware(**cfg.daq_obj_kwds)
# ni.configure(cfg.get_daq_cycle_time(), cfg.daq_ao_names_to_channels, live = True)
# _, voltages_t = generate_waveforms(cfg, 488)
# ni.assign_waveforms(voltages_t)

# Setup Camera
runtime = calliphlox.Runtime()
logging.info('Initializing Camera')
dm = runtime.device_manager()
p = runtime.get_configuration()

cameras = [
    d.name
    for d in dm.devices()
    if (d.kind == DeviceKind.Camera) and ("C15440" in d.name)
]
p.video[0].camera.identifier = dm.select(DeviceKind.Camera, cameras[0])
p.video[0].camera.settings.exposure_time_us = 0
#p.video[0].camera.settings.line_interval_us = 25
p.video[0].storage.settings.filename = fr"C:\Acquire Test_0_Zarr.zarr"  # Need to put something at the beginning but can change later
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

# Setup log
FILEPATH = r"C:\Acquire Test\1_zarr_testing_02_2_23.log"
logging.getLogger().setLevel(logging.DEBUG)
# Create a file handler.
log_handler = logging.FileHandler(FILEPATH, 'w')
log_handler.setLevel(logging.DEBUG)
logging.getLogger().addHandler(log_handler)

frequencys = [30]#,40,50,60,70,80,90]
exposure_times = [20]#,10,15,20,25]
zarr_types = ["Tiff"] #"Zarr", "ZarrBlosc1ZstdByteShuffle"]
for freqs in frequencys:

    ni.livestream_frequency_hz = freqs
    logging.info(f"Setting frequency to {ni.livestream_frequency_hz}")
    ni.configure(cfg.get_daq_cycle_time(), cfg.daq_ao_names_to_channels, live=True)
    _, voltages_t = generate_waveforms(cfg, 488)
    ni.assign_waveforms(voltages_t)

    for time in exposure_times:

        logging.info(f"Setting exposure time to {time}")
        logging.info(f"Setting linerate to 0")
        p.video[0].camera.settings.exposure_time_us = time
        #p.video[0].camera.settings.line_interval_us = 25
        logging.info(f"Exposure time set to {p.video[0].camera.settings.exposure_time_us}")
        logging.info(f"Linerate set to {p.video[0].camera.settings.line_interval_us}")

        for type in zarr_types:
            logging.info(f"Setting file saving to {type}")
            filename = fr"C:\Acquire Test_exposure_{time}_frequency_{freqs}_{type}"
            if type == 'Tiff':
                try:
                    os.mkdir(filename)
                except:
                    pass
                filepath = fr"{filename}\{filename[3:-1]}.tiff"
            else:
                filepath = f"{filename}.zarr"


            p.video[0].storage.settings.filename = filepath
            p.video[0].storage.identifier = dm.select(DeviceKind.Storage, type)
            runtime.set_configuration(p)

            logging.info("Starting cameras.")
            print(p.video[0])
            runtime.start()

            frames_collected = 0
            logging.info("Starting NI card")
            ni.start()
            ni.ci_task.start()

            timeout = now() + (1000/freqs)  #How long it will take to collect 1000 frames
            logging.info(f"Run set for {1000/freqs} seconds")

            #while timeout > now():
            counts = 0
            while counts < 1000:
                counts = ni.ci_task.read()
                if a := runtime.get_available_data(0):
                    packet = a.get_frame_count()
                    f = next(a.frames())
                    im = f.data().squeeze().copy()

                    f = None  # <-- fails to get the last frames if this is held?
                    a = None  # <-- fails to get the last frames if this is held?

                    frames_collected += packet
                    #print(f'frames collected: {frames_collected}')


            logging.info(f'Run ended. {frames_collected} out of 1000')

            logging.info("Stopping NI card")
            ni.stop()
            ni.ci_task.stop()
            # ni.close()
            # ni.ci_task.close()

            runtime.abort()
            logging.info('Camera stop')

            logging.info(f'Completed run exposure time: {time} and file storage: {type} at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")} ')


runtime.abort()
ni.close()
print('Finished')

