import nidaqmx
import numpy as np
from nidaqmx.constants import AcquisitionType as AcqType
from nidaqmx.constants import TaskMode, FrequencyUnits, Level
from nidaqmx.constants import Edge, Slope
from numpy import ndarray
from dispim.dispim_config import DispimConfig
from dispim.devices.ni import WaveformHardware
from dispim.compute_waveforms import generate_waveforms
import time

cfg = DispimConfig(r'C:\Users\Administrator\Projects\dispim-control\examples\config.toml')
ni = WaveformHardware(**cfg.daq_obj_kwds)
ni.configure(cfg.get_daq_cycle_time(), cfg.daq_ao_names_to_channels, live = True)
_, voltages_t = generate_waveforms(cfg, 488)
ni.assign_waveforms(voltages_t)
ni.start()
time.sleep(5)

# counter_task = nidaqmx.Task("counter_task")
# counter_task.co_channels.add_co_pulse_chan_freq('/Dev2/ctr0',
#                                                                  units=FrequencyUnits.HZ,
#                                                                  idle_state=Level.LOW, initial_delay=0.0,
#                                                                  freq=15,  # change 15 - 30 Hz, change to config value
#                                                                  duty_cycle=0.5)
# counter_task.ci_count_edges_term = '/Dev2/PFI3'
# counter_task.timing.cfg_implicit_timing(sample_mode=AcqType.CONTINUOUS)
# counter_task.start()
# time.sleep(10)