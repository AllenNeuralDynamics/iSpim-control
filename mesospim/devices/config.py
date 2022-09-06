"""Read config file and compute derived parameters."""

import toml
from math import ceil

# CONSTANTS:
# for sanity checks. Do not change these.
MIN_LINE_INTERVAL = 10e-6  # seconds. min time a row is turned on.
MAX_GALVO_FREQ_HZ = 1000
SENSOR_ROW_COUNT = 2048
MAX_ADJ_GALVO_VOLTAGE = 10 # max differential voltage abs(a-b) that the galvo can withstand.


# Retrieve the config file.
cfg_file = None
with open('config.toml', 'r') as toml_file:
    cfg_file = toml.load(toml_file)


# Tunable Parameters
slit_width = cfg_file['experiment_params']['slit_width']  # 274 nominal.
row_interval = cfg_file['dcam_params']['row_interval']  # Time between moving the slit by one row.
                                                        # aka "line interval."

# Helper value
row_exposure_time = row_interval * slit_width  # total time a row gets exposed to laser.
total_exposure_time = SENSOR_ROW_COUNT*row_interval + row_exposure_time  # total time spent
                                                                         # where **any** row is being exposed.

# Derived or Fixed Parameters
galvo_freq = 1/row_exposure_time  # To set galvo freq, we must sweep at least once in the time a row is turned on.

etl_settling_time = cfg_file['waveform_specs']['etl']['settling_time']
# recovery_time is extra time for all devices to return to their start positions
# Make recovery_time
# 1. an integer multiple of galvo periods.
# 2. at least as long as the etl_settling_time
recovery_time = ceil(total_exposure_time*galvo_freq)/galvo_freq - total_exposure_time  # remainder.
while recovery_time < etl_settling_time:
    recovery_time += row_exposure_time
waveform_cycle_time = total_exposure_time + recovery_time # cycle time of waveforms not including delays.

# Waveform Specs. Pull most from the config. Derive the remaining ones.
etl_delay = cfg_file['waveform_specs']['etl']['delay']
start_of_frame_delay = cfg_file['waveform_specs']['start_of_frame_delay']
waveform_delay = etl_delay + start_of_frame_delay # delay for all waveforms (except etl)
daq_cycle_time = waveform_cycle_time + waveform_delay
galvo_min_voltage = cfg_file['waveform_specs']['galvo']['min_voltage']
galvo_max_voltage = cfg_file['waveform_specs']['galvo']['max_voltage']

etl_min_voltage = cfg_file['waveform_specs']['etl']['min_voltage']
etl_max_voltage = cfg_file['waveform_specs']['etl']['max_voltage']
etl_duty = 1.0 # pure sawtooth
#etl_duty = total_exposure_time / waveform_cycle_time

laser_enable_min_voltage = cfg_file['waveform_specs']['laser_enable']['min_voltage']
laser_enable_max_voltage = cfg_file['waveform_specs']['laser_enable']['max_voltage']
laser_enable_duty = total_exposure_time / waveform_cycle_time

adj_galvo_setpoint = cfg_file['waveform_specs']['adjustment_galvo']['voltage_setpoint']

# Daq Hardware Specs
daq_used_channels = len(cfg_file['daq_params']['ao_names_to_channels']) + 1 # add 1 for digital out.
daq_update_freq = cfg_file['daq_params']['rate']  # 400 KHz nominal.
daq_num_samples = round(daq_update_freq*daq_cycle_time)

# Sanity Checks
assert row_interval >= MIN_LINE_INTERVAL, \
    f"Error: row interval ({row_interval*1e6:.3f} [us]) too small. Minimum is {MIN_LINE_INTERVAL*1e6:.3f} [us]."
assert galvo_freq <= MAX_GALVO_FREQ_HZ or (galvo_min_voltage == galvo_max_voltage), \
    f"Error: desired galvo frequency ({galvo_freq} [Hz]) too fast. Max is ({MAX_GALVO_FREQ_HZ} [Hz])."
assert abs(adj_galvo_setpoint * 2) < MAX_ADJ_GALVO_VOLTAGE, \
    f"Error: adjustment galvo setpoint voltage is out of range."


# EXPERIMENT PARAMETERS loaded from file.
path = cfg_file['experiment_params']['path']
filename = cfg_file['experiment_params']['filename']  # specimen names
#nFrames = cfg_file['experiment_params']['frames_to_grab']  # number of frames to grab

tile_size_um = cfg_file['experiment_params']['tile_size_um']
iso_voxel_size_um = tile_size_um / SENSOR_ROW_COUNT # size of one voxel cube.
percent_image_overlap = cfg_file['experiment_params']['percent_image_overlap']
grid_step_um = (1 - percent_image_overlap/100.0)*tile_size_um

volume_x_um = cfg_file['experiment_params']['volume_x_um']
volume_y_um = cfg_file['experiment_params']['volume_y_um']
volume_z_um = cfg_file['experiment_params']['volume_z_um']

x_steps = ceil((volume_x_um - tile_size_um)/grid_step_um)
y_steps = ceil((volume_y_um - tile_size_um)/grid_step_um)
z_steps = ceil((volume_z_um - iso_voxel_size_um)/iso_voxel_size_um)

x_tiles, y_tiles, z_tiles = (1+x_steps, 1+y_steps, 1+z_steps)
total_tiles = x_tiles*y_tiles*z_tiles

actual_volume_x_um = tile_size_um + x_steps*grid_step_um
actual_volume_y_um = tile_size_um + y_steps*grid_step_um
actual_volume_z_um = iso_voxel_size_um*(1 + z_steps)

row_interval = cfg_file['dcam_params']['row_interval']
slit_width = cfg_file['experiment_params']['slit_width']

# DEVICE PARAMETERS loaded from file.
scan_direction = cfg_file['dcam_params']['scan_direction']
etl_params = cfg_file['etl_params']
daq_params = cfg_file['daq_params']
daq_params['periodTime'] = daq_cycle_time  # add computed parameter.


# Print some stats:
print("--Config Stats--")
print(f"  Galvo frequency (computed): {galvo_freq:.3f} [Hz]")
print(f"  Total \"sensor-active\" time (computed): {total_exposure_time*1e3:.3f} [ms]")
print(f"  DAQ cycle time (computed): {daq_cycle_time*1e3:.3f} [ms]")
print()
print("Volume Capture Stat:")
print(f"  grid step: {grid_step_um} [um]")
print(f"  percent overlap: {percent_image_overlap:.1f}%")
print(f"  Desired dimensions: {volume_x_um:.1f}[um] x {volume_y_um:.1f}[um] x {volume_z_um:.1f}[um]")
print(f"  Actual dimensions:  {actual_volume_x_um:.1f}[um] x {actual_volume_y_um:.1f}[um] x {actual_volume_z_um:.1f}[um]")
print(f"  Actual dimensions:  {x_tiles}[tiles] x {y_tiles}[tiles] x {z_tiles}[tiles]")
print(f"  Total number of tiles: {total_tiles}")
print(f"  voxel size: {iso_voxel_size_um:.3f}[um]")
print()
print("GUI settings for debugging:")
print(f"  Line Interval: {row_interval*1e6:.3f} [us]")
print(f"  Exposure Time (computed): {row_exposure_time*1e3:.3f} [ms]")
print()
