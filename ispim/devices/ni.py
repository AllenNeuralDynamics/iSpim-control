"""
NIDAQMX setup for timing/triggering using PCIe-6738 card.
"""

import nidaqmx
import numpy as np
from nidaqmx.constants import AcquisitionType as AcqType
from nidaqmx.constants import TaskMode, FrequencyUnits, Level
from nidaqmx.constants import Edge, Slope
from numpy import ndarray
import logging
from time import sleep


class WaveformHardware:
    def __init__(self, dev_name, input_trigger_name, count_trigger_name,
                 ao_counter_trigger_name, update_frequency_hz, livestream_frequency_hz):

        self.log = logging.getLogger(__name__)
        self.dev_name = dev_name  # NI card address, i.e. Dev2
        self.dev = nidaqmx.system.device.Device(self.dev_name)
        self.log.warning('Resetting NIDAQ')
        self.dev.reset_device()
        self.input_trigger_name = input_trigger_name.lstrip('/')  # NI card output trigger port, i.e. PFI00
        self.counter_trigger_name = count_trigger_name.lstrip('/')  # PFI counter task port
        self.ao_counter_trigger_name = ao_counter_trigger_name.lstrip('/')  # PFI port that triggers ao lines (connected to PFI counter task port)
        self.update_freq = update_frequency_hz  # in [Hz]
        self.livestream_frequency_hz = livestream_frequency_hz
        self.ao_task = None
        self.counter_task = None
        self.live = None

    def configure(self, period_time: float, ao_names_to_channels: dict, channel_num : int = 1,
                  live: bool = False):
        """Configure the daq with tasks."""
        # Close any existing tasks if we are reconfiguring.
        if self.counter_task or self.ao_task:
            if not self.live:
                self.wait_until_done()
            self.stop()
            self.close()
        self.live = live
        sample_count = round(self.update_freq * period_time)    # sample_count for single channel. All the same sample per channel
        # Create AO task and initialize the required channels
        self.ao_task = nidaqmx.Task("analog_output_task")
        for channel_name, channel_index in ao_names_to_channels.items():
            physical_name = f"/{self.dev_name}/ao{channel_index}"
            self.log.debug(f"Setting up ao channel {channel_name} "
                            f"on {physical_name}")
            self.ao_task.ao_channels.add_ao_voltage_chan(physical_name)
        self.ao_task.timing.cfg_samp_clk_timing(
            rate=self.update_freq,
            active_edge=Edge.RISING,
            sample_mode=AcqType.FINITE,
            samps_per_chan=sample_count*channel_num)
        # Takes into account the number of channels in ao task
        self.ao_task.triggers.start_trigger.retriggerable = True

        if live:
            self.counter_task = nidaqmx.Task("counter_task")
            co_channel = self.counter_task.co_channels.add_co_pulse_chan_freq(f"/{self.dev_name}/ctr0",
                                                            units=FrequencyUnits.HZ,
                                                            idle_state=Level.LOW, initial_delay=0.0,
                                                            freq= self.livestream_frequency_hz,  # change 15 - 30 Hz, change to config value
                                                            duty_cycle=0.5)
            co_channel.co_pulse_term = f"/{self.dev_name}/{self.counter_trigger_name}"
            self.counter_task.timing.cfg_implicit_timing(sample_mode=AcqType.CONTINUOUS)

            self.ao_task.triggers.start_trigger.cfg_dig_edge_start_trig(trigger_source=f"/{self.dev_name}/{self.ao_counter_trigger_name}",
                                                                        # if in live mode PFI3 trigger_edge = Slope.RISING)
                                                                        trigger_edge=Slope.RISING)


        else:

            self.ao_task.triggers.start_trigger.cfg_dig_edge_start_trig(
                trigger_source=f"/{self.dev_name}/{self.input_trigger_name}",
                trigger_edge=Slope.RISING)
            
            self.counter_task = nidaqmx.Task("counter_task")
            self.counter_loop = self.counter_task.ci_channels.add_ci_count_edges_chan(f'/{self.dev_name}/ctr0',
                                                                                      edge=nidaqmx.constants.Edge.RISING)
            self.counter_loop.ci_count_edges_term = f"/{self.dev_name}/{self.input_trigger_name}"

            # NOT SURE IF WE NEED THIS?
            # "Commit" if we're not looping. Apparently, this has less overhead.
            # https://forums.ni.com/t5/LabVIEW/Deleting-channels-from-task-reconfiguring-task/m-p/1544490/highlight/true#M571637

            self.ao_task.out_stream.output_buf_size = sample_count*channel_num  # Sets buffer to length of voltages
            self.ao_task.control(TaskMode.TASK_COMMIT)

    def assign_waveforms(self, voltages_t, scout_mode: bool = False):
        """Write analog and digital waveforms to device.
        Order is driven by the TOML config file.
        """
        # Confirm digital signal waveform is a numpy array because we must
        # ultimately write the digital waveform to the device as bools.

        assert type(voltages_t) == ndarray, \
            "Error: voltages_t digital signal waveform must be a numpy ndarray."
        # Write analog voltages.
        if scout_mode:
            self.rereserve_buffer(len(voltages_t[0]))
        self.ao_task.write(voltages_t, auto_start=False)  # arrays of floats


    def start(self):
        """start tasks."""
        # Start ao_task and counter task
        self.ao_task.start()
        self.counter_task.start()

    def playback_finished(self):
        """True if the device is busy playing waveforms. False otherwise."""
        # Check if ao task is finished
        return self.ao_task.is_task_done()

    def wait_until_done(self, timeout=1.0):
        # Check if ao task is finished
        return self.ao_task.wait_until_done(timeout)

    def stop(self, wait = 1):
        """Stop the tasks"""
        # For tasks in continuous mode, we need to write all zeros before
        # closing to ensure the waveforms exit in known state.
        # We need to write an array as the same size as what the task was
        # provisioned with.
        # if self.live:
        #     ao_data = np.zeros((len(self.ao_task.ao_channels), self.ao_task.out_stream.output_buf_size))
        #     self.ao_task.write(ao_data)
        #TODO: Why is this not working?
        self.log.debug("Issuing a task stop.")
        self.counter_task.stop()
        self.counter_task.wait_until_done(1)
        sleep(wait)         # Sleep so ao task can finish
        self.ao_task.stop()

    def rereserve_buffer(self, buf_len):
        """If tasks are already configured, the buffer needs to be cleared and rereserved to work"""

        if self.ao_task is not None:
            self.ao_task.control(TaskMode.TASK_UNRESERVE)  # Unreserve buffer
            self.ao_task.out_stream.output_buf_size = buf_len  # Sets buffer to length of voltages
            self.ao_task.control(TaskMode.TASK_COMMIT)


    def restart(self):
        # Restart ao task
        self.ao_task.stop()
        self.ao_task.start()

    def close(self):
        """Terminate all started tasks."""
        if self.counter_task is not None:
            self.log.debug("closing tasks!")
            self.counter_task.close()
            self.counter_task = None
        if self.ao_task is not None:
            self.ao_task.close()
            self.ao_task = None
