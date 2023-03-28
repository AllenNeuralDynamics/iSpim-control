import nidaqmx
import numpy as np
from nidaqmx.constants import AcquisitionType
from nidaqmx.constants import TaskMode, FrequencyUnits, Level
from nidaqmx.constants import Edge, Slope
from time import sleep
from calliphlox import DeviceKind, Trigger, SampleType, TriggerEvent, SignalIOKind, TriggerEdge, Direction
import calliphlox
import logging
import os
from datetime import datetime

class NI:

    def configure_tasks(self, frequency):

        # Setting up counter task to trigger ao task
        self.co_task = nidaqmx.Task("co_task")
        co_channel = self.co_task.co_channels.add_co_pulse_chan_freq('/Dev2/ctr0',units=FrequencyUnits.HZ,
                                                                    idle_state=Level.LOW, initial_delay=0.0,
                                                                    freq= frequency,  # change 15 - 30 Hz, change to config value
                                                                    duty_cycle=0.5)

        co_channel.co_pulse_term = '/Dev2/PFI3'
        self.co_task.timing.cfg_implicit_timing(sample_mode=AcquisitionType.CONTINUOUS)

        # Setting up counter task to count pulses

        self.ci_task = nidaqmx.Task("ci_task")
        ci_channel = self.ci_task.ci_channels.add_ci_count_edges_chan('/Dev2/ctr1', edge = nidaqmx.constants.Edge.RISING)
        ci_channel.ci_count_edges_term = '/Dev2/PFI3'

        # Setting up ao task to trigger camera
        sample_count = round(400000.0 * (0.006+0.004+0.025))    # Update frequency plus delay, exposure, and rest time
        self.ao_task = nidaqmx.Task("analog_output_task")
        self.ao_task.ao_channels.add_ao_voltage_chan('/Dev2/ao6')
        self.ao_task.timing.cfg_samp_clk_timing(
            rate=400000.0,
            active_edge=Edge.RISING,
            sample_mode=AcquisitionType.FINITE,
            samps_per_chan=sample_count)
        self.ao_task.triggers.start_trigger.retriggerable = True
        self.ao_task.triggers.start_trigger.cfg_dig_edge_start_trig(
            trigger_source='/Dev2/PFI2',
            trigger_edge=Slope.RISING)


    def write_waveforms(self):

        voltages_t = np.zeros(round(400000.0 * (0.006+0.004+0.025)))
        voltages_t[round(400000.0 * 0.006) + round(400000.0 * 0.001):
                   round(400000.0 * 0.006) + round(400000.0 * 0.001) + round(400000.0 * 0.025)] = 5.0

        self.ao_task.write(voltages_t, auto_start=False)

    def start(self):

        self.co_task.start()
        self.ci_task.start()
        self.ao_task.start()

    def stop(self):

        self.co_task.stop()
        self.ci_task.stop()
        self.ao_task.stop()

    def close(self):
        self.co_task.close()
        self.ci_task.close()
        self.ao_task.close()

class FrameGrabber:

    def initialize_camera(self, runtime):
        dm = runtime.device_manager()
        p = runtime.get_configuration()

        cameras = [
            d.name
            for d in dm.devices()
            if (d.kind == DeviceKind.Camera) and ("C15440" in d.name)
        ]

        return p, cameras

    def setup_camera(self, runtime, p, cameras):
        dm = runtime.device_manager()
        p.video[0].camera.identifier = dm.select(DeviceKind.Camera, cameras[0])
        p.video[0].camera.settings.exposure_time_us = 4340.28
        p.video[0].camera.settings.line_interval_us = 10.85
        p.video[0].storage.identifier = dm.select(DeviceKind.Storage, 'Trash')
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

        runtime.set_configuration(p)
        return runtime, p, dm

if __name__ == "__main__":

    # Initializing NI Class
    ni = NI()
    ni.configure_tasks(30)
    ni.write_waveforms()

    # ni.configure_tasks(30)
    # ni.write_waveforms()
    # ni.start()
    # counts = 0
    # while counts < 1000:
    #     counts = ni.ci_task.read()
    #     print(counts)
    #     sleep(0.01)
    # ni.stop()
    # ni.close()

    #Initializing framegrabber
    framegrabber = FrameGrabber()
    runtime = calliphlox.Runtime()

    logging.info('Initializing Camera')
    p, cameras = framegrabber.initialize_camera(runtime)
    runtime, p, dm = framegrabber.setup_camera(runtime, p, cameras)

    # Initializing logging
    FILEPATH = r"C:\Acquire Test\single_camera_test_03_15_23.log"
    logging.getLogger().setLevel(logging.DEBUG)
    # Create a file handler.
    log_handler = logging.FileHandler(FILEPATH, 'w')
    log_handler.setLevel(logging.DEBUG)
    logging.getLogger().addHandler(log_handler)


    while True:

        runtime.start()

        frames_collected = 0
        logging.info("Starting NI card")
        ni.start()
        counts = 0
        while counts < 20000:
            counts = ni.ci_task.read()
            if a := runtime.get_available_data(0):
                packet = a.get_frame_count()
                f = next(a.frames())
                im = f.data().squeeze().copy()

                f = None  # <-- fails to get the last frames if this is held?
                a = None  # <-- fails to get the last frames if this is held?

                frames_collected += packet
                print(f'frames collected: {frames_collected}')

        logging.info(f'Run ended. \033[1m{frames_collected}\033[0m out of 20000')

        logging.info("Stopping NI card")
        ni.stop()

        runtime.abort()
        logging.info('Camera stop')

        logging.info(
            f'Completed run at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")} ')

    ni.close()

    runtime.abort()
    ni.close()
    print('Finished')


