"""
NIDAQMX setup for timing/triggering using PCIe-6738 card.

"""

import nidaqmx
import numpy as np
from nidaqmx.constants import AcquisitionType as AcqType
from nidaqmx.constants import TaskMode
from nidaqmx.constants import Edge, Slope
from numpy import ndarray
import logging


class WaveformHardware:

    def __init__(self, dev_name, output_trigger_name, update_frequency):

        self.dev_name = dev_name # NI card address, i.e. Dev2
        self.output_trigger_name = output_trigger_name.lstrip('/')  # NI card output trigger port, i.e. PFI00
        self.update_freq = update_frequency  # in [Hz]
        self.ao_task = None
        self.do_task = None
        self.log = logging.getLogger(__name__)
        self.live = None

    def configure(self, period_time: float, ao_names_to_channels: dict,
                  live: bool = False):
        """Configure the daq with tasks."""
        # Close any existing tasks if we are reconfiguring.
        if self.do_task or self.ao_task:
            if not self.live:
                self.wait_until_done()
            self.stop()
            self.close()
        self.live = live
        sample_count = round(self.update_freq*period_time)
        # Create a Daq DO and AO tasks and initialize the required channels.
        # Configure digital output task.
        self.do_task = nidaqmx.Task("digital_output_task")
        do_physical_name = f"/{self.dev_name}/{self.output_trigger_name}"
        self.log.debug(f"Setting up camera external trigger "
                       f"on {do_physical_name}")
        self.do_task.do_channels.add_do_chan(do_physical_name)
        self.do_task.timing.cfg_samp_clk_timing(
            rate=self.update_freq,
            active_edge=Edge.RISING,
            sample_mode=AcqType.CONTINUOUS if live else AcqType.FINITE,
            samps_per_chan=sample_count)
        # Digital output task will trigger when we call start().

        # Configure analog output task. Set triggering to the start of digital output task. See:
        # https://www.ni.com/docs/en-US/bundle/ni-daqmx-21.3-help/page/mxcncpts/syncstarttrigger.html
        self.ao_task = nidaqmx.Task("analog_output_task")
        for channel_name, channel_index in ao_names_to_channels.items():
            physical_name = f"/{self.dev_name}/ao{channel_index}"
            self.log.debug(f"Setting up ao channel {channel_name} "
                           f"on {physical_name}")
            self.ao_task.ao_channels.add_ao_voltage_chan(physical_name)
        self.ao_task.timing.cfg_samp_clk_timing(
            rate=self.update_freq,
            active_edge=Edge.RISING,
            sample_mode=AcqType.CONTINUOUS if live else AcqType.FINITE,
            samps_per_chan=sample_count)
        self.ao_task.triggers.start_trigger.retriggerable = not live
        # Configure ao_task to start each time the do_task starts.
        self.ao_task.triggers.start_trigger.cfg_dig_edge_start_trig(
            trigger_source=f"/{self.dev_name}/do/StartTrigger",
            trigger_edge=Slope.RISING)
        # "Commit" if we're not looping. Apparently, this has less overhead.
        # https://forums.ni.com/t5/LabVIEW/Deleting-channels-from-task-reconfiguring-task/m-p/1544490/highlight/true#M571637
        self.ao_task.out_stream.output_buf_size = sample_count
        self.ao_task.control(TaskMode.TASK_COMMIT)
        self.do_task.out_stream.output_buf_size = sample_count
        self.do_task.control(TaskMode.TASK_COMMIT)

    def assign_waveforms(self, voltages_t):
        """Write analog and digital waveforms to device.
        Order is driven by the TOML config file.
        """
        # Confirm digital signal waveform is a numpy array because we must
        # ultimately write the digital waveform to the device as bools.

        assert type(voltages_t[-1]) == ndarray, \
            "Error: voltages_t digital signal waveform must be a numpy ndarray."
        # Write analog voltages.
        self.ao_task.write(voltages_t[:-1], auto_start=False)  # arrays of floats
        # Write the digital output trigger signal.
        self.do_task.write(voltages_t[-1].astype(bool), auto_start=False)  # an array of bools.

    def start(self):
        """start tasks."""
        # Start ao_task first since do_task's start will act as trigger source for the ao task.
        self.ao_task.start()
        self.do_task.start()

    def playback_finished(self):
        """True if the device is busy playing waveforms. False otherwise."""
        # Only check do_task since ao_task is subordinate to it.
        return self.do_task.is_task_done()

    def wait_until_done(self, timeout=1.0):
        # Only check do_task since ao_task is subordinate to it.
        return self.do_task.wait_until_done(timeout)

    def stop(self):
        """Stop the tasks"""
        # For tasks in continuous mode, we need to write all zeros before
        # closing to ensure the waveforms exit in known state.
        # We need to write an array as the same size as what the task was
        # provisioned with.
        if self.live:
            ao_data = np.zeros((len(self.ao_task.ao_channels), self.ao_task.out_stream.output_buf_size))
            self.ao_task.write(ao_data)
            self.do_task.write([False]*self.do_task.out_stream.output_buf_size,
                               auto_start=True)
        # Stop the DO task first since it is the trigger and we
        # don't want to kill waveform playback prematurely.
        self.log.debug("Issuing a task stop.")
        self.do_task.stop()
        self.ao_task.stop()

    def restart(self):
        # ao_task is configured to start when the do task starts.
        self.do_task.stop()
        self.do_task.start()

    def close(self):
        """Terminate all started tasks."""
        if self.do_task is not None:
            self.log.info("closing tasks!")
            self.do_task.close()
            self.do_task = None
        if self.ao_task is not None:
            self.ao_task.close()
            self.ao_task = None
