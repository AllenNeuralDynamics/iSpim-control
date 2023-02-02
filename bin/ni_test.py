import nidaqmx
import numpy as np
from nidaqmx.constants import AcquisitionType 
from nidaqmx.constants import TaskMode, FrequencyUnits, Level
from nidaqmx.constants import Edge, Slope
from numpy import ndarray
from dispim.dispim_config import DispimConfig
from dispim.devices.ni import WaveformHardware
from dispim.compute_waveforms import generate_waveforms
import time
#
# cfg = DispimConfig(r'C:\Users\Administrator\Projects\dispim-control\examples\config.toml')
# ni = WaveformHardware(**cfg.daq_obj_kwds)
# ni.configure(cfg.get_daq_cycle_time(), cfg.daq_ao_names_to_channels, live = True)
# _, voltages_t = generate_waveforms(cfg, 488)
# ni.assign_waveforms(voltages_t)
# ni.start()
# time.sleep(5)

co_task = nidaqmx.Task("co_task")
co_channel = co_task.co_channels.add_co_pulse_chan_freq('/Dev2/ctr0',units=FrequencyUnits.HZ,
                                                            idle_state=Level.LOW, initial_delay=0.0,
                                                            freq= 15,  # change 15 - 30 Hz, change to config value
                                                            duty_cycle=0.5)

co_channel.co_pulse_term = '/Dev2/PFI3'
co_task.timing.cfg_implicit_timing(sample_mode=AcquisitionType.CONTINUOUS)

ci_task = nidaqmx.Task("ci_task")
ci_channel = ci_task.ci_channels.add_ci_count_edges_chan('/Dev2/ctr1', edge = nidaqmx.constants.Edge.RISING)
ci_channel.ci_count_edges_term = '/Dev2/PFI3'

co_task.start()
ci_task.start()
counts = 0
while counts < 10000:
    counts = ci_task.read()
    print(counts)
    time.sleep(0.01)
