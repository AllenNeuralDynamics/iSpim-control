#!/usr/bin/env python3
"""Generate waveforms from config params. If standalone, save a graph only."""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from scipy.signal import square, sawtooth
from .mesospim_config import MesospimConfig

# Notes:
# https://github.com/mesoSPIM/mesoSPIM-hardware-documentation/wiki/mesoSPIM_waveforms


# TODO: cfg should be able to lookup config params sensibly (like with a string for a key)
def generate_waveforms(cfg: MesospimConfig, active_wavelength: int):
    """return a numpy nd array with the correct waveforms.

    The DAQ outputs all waveforms to control: the etl, the galvo,
    laser enable(s), and the start trigger.
    The camera will begin capturing upon reading the start trigger.

    Waveforms are planned to account for any non-negligible phase lag.
    Here, the ETL is laggiest. ETL signal will start first such that the true
    ETL output occurs in sync with the remaining signals.

    Since multiple lasers are present, a waveform is generated for all
    waveforms, but only the active ones will have nonzero values.

    :param cfg: Mesospim configuration object.
    :param active_wavelength: laser wavelength that will be turned on.
    """
    # Create wavelength-dependent constants
    active_laser_specs = cfg.laser_specs[str(active_wavelength)]['etl']
    etl_min_voltage = active_laser_specs['min_voltage']
    etl_max_voltage = active_laser_specs['max_voltage']
    waveform_cycle_time = cfg.get_waveform_cycle_time(active_wavelength)
    etl_delay = cfg.get_etl_delay(active_wavelength)
    daq_cycle_time = cfg.get_daq_cycle_time(active_wavelength)
    daq_num_samples = cfg.get_daq_num_samples(active_wavelength)
    laser_duty_cycle = cfg.get_laser_duty_cycle(active_wavelength)

    # Create table that holds an entire cycle's worth of pts.
    # External trigger signal is the last value. Aside from that,
    # the signal order in voltages_t doesn't matter provide that it is
    # consistent with the signal order that the WaveformGenerator creates
    # tasks. We can ensure this by iterating through the same data structure
    # in the cfg as we do when we create tasks with the WaveformGenerator.
    voltages_t = np.zeros((cfg.daq_used_channels, daq_num_samples))
    t = np.linspace(0, daq_cycle_time, daq_num_samples, endpoint=False)

    # Delay for all signals except the external trigger.
    sof_shift_count = round(cfg.daq_update_freq * cfg.start_of_frame_delay)
    # Delay for all signals except the ETL.
    etl_shift_count = round(cfg.daq_update_freq * etl_delay)
    # Iterate through all laser channels and create the corresponding waveform.
    for index, ao_name in enumerate(cfg.daq_ao_names_to_channels.keys()):
        # Generate Laser Signal. Signal should be:
        # a pulse for the active laser but a flatline for inactive lasers.
        if ao_name in [str(nm) for nm in cfg.laser_wavelengths]:
            specs = cfg.laser_specs[ao_name]
            disable_voltage = specs['disable_voltage']
            enable_voltage = specs['disable_voltage']
            # Only enable the active wavelengths.
            if ao_name == str(active_wavelength):
                enable_voltage = specs['enable_voltage']
            # Generate Laser Signal analog time series.
            laser_amplitude = (enable_voltage - disable_voltage)/2.0
            laser_dc_offset = disable_voltage + laser_amplitude
            voltages_t[index] = laser_amplitude * \
                                square(2*pi / waveform_cycle_time*t,
                                       duty=laser_duty_cycle) + laser_dc_offset
            # Apply start of frame delay to all waveforms (except external trigger)
            voltages_t[index] = np.roll(voltages_t[index], sof_shift_count)
            # Stuff the trough value at the beginning.
            np.put(voltages_t[index], range(sof_shift_count),
                   laser_dc_offset - laser_amplitude)
            # Apply ETL shift
            voltages_t[index] = np.roll(voltages_t[index], etl_shift_count)
            # Stuff the trough value at the beginning.
            np.put(voltages_t[index], range(etl_shift_count),
                   laser_dc_offset - laser_amplitude)
        # Generate ETL signal, a pure sawtooth followed by delay time.
        elif ao_name.startswith("etl"):
            # ETL settings are specific to the active laser.
            etl_amplitude = (etl_max_voltage - etl_min_voltage) / 2.0
            etl_offset = etl_min_voltage + etl_amplitude
            sawtooth_period = cfg.total_exposure_time + etl_delay
            voltages_t[index] = \
                etl_amplitude * sawtooth(2*pi/sawtooth_period*t, width=1.0) + \
                etl_offset
            # Flatten out end of sawtooth wave by stuffing it with min value.
            end_of_sawtooth_sample = round(sawtooth_period*cfg.daq_update_freq)
            np.put(voltages_t[index],
                   range(end_of_sawtooth_sample, daq_num_samples),
                   etl_offset - etl_amplitude)
            # Apply start of frame delay to all waveforms (except external trigger)
            voltages_t[index] = np.roll(voltages_t[index], sof_shift_count)
            # Stuff the trough value at the beginning.
            np.put(voltages_t[index], range(sof_shift_count),
                   etl_offset - etl_amplitude)
        # Generate differential voltages for the adjustment galvo.
        # These signals are fixed throughout and don't need shifting.
        elif ao_name.startswith("galvo") and ao_name.endswith("plus"):
            voltages_t[index] = 0.5 * cfg.adj_galvo_setpoint * np.ones(len(t))
        elif ao_name.startswith("galvo") and ao_name.endswith("minus"):
            voltages_t[index] = -0.5 * cfg.adj_galvo_setpoint * np.ones(len(t))
        else:
            raise RuntimeError(f"{ao_name} does not have any plotting criteria.")

    # Generate External Trigger signal.
    np.put(voltages_t[-1], range(round(10e-6 * cfg.daq_update_freq)), 1)
    # Apply etl delay.
    voltages_t[-1] = np.roll(voltages_t[-1], etl_shift_count)
    np.put(voltages_t[-1], range(etl_shift_count), 0)

    return t, voltages_t


def plot_waveforms_to_pdf(cfg: MesospimConfig, t: list, voltages_t: list,
                          active_wavelength: int, filename: str = "plot.pdf"):
    """Write a pdf plot output of the waveforms."""
    # Plot the data for sanity checking.
    fig = plt.figure(figsize=(20, 7))
    ax = fig.gca()
    # first plot: the whole thing.
    fig.suptitle("One Image Capture Sequence.")
    for index, ao_name in enumerate(cfg.daq_ao_names_to_channels.keys()):
        ax.plot(t, voltages_t[index], label=ao_name)
    # Last time series is the Camera Trigger and is not an analog output.
    ax.plot(t, voltages_t[-1], label="camera trigger")
    # Plot desired etl signal
    #etl_shift_count = round(cfg.daq_update_freq * cfg.etl_delay)
    #axes[0].plot(t, np.roll(etl_t, etl_shift_count), label="desired liquid lens")

    # Include exposure time stats on plot.
    last_row_start_time = cfg.row_interval*cfg.sensor_row_count + \
                          cfg.get_waveform_delay(active_wavelength)
    last_row_stop_time = cfg.total_exposure_time + \
                         cfg.get_waveform_delay(active_wavelength)
    ax.axvline(last_row_start_time, ls='--', color='cyan',
                    label=f"last row start: {last_row_start_time*1e3:.2f}[ms]")
    ax.axvline(last_row_stop_time, ls='--', color='magenta',
                    label=f"last row fin: {last_row_stop_time * 1e3:.2f}[ms]")

    ax.set_xlabel("time [s]")
    ax.set_ylabel("amplitude [V]")
    ax.legend(loc="upper right")#loc="center left")

    try:
        fig.savefig(filename)
    except OSError as e:
        print("Error: cannot save figure. Another program may be using it.")
        raise e


if __name__ == "__main__":
    from argparse import ArgumentParser
    from .config_base import MesospimConfig
    # Argparse for a config file.
    parser = ArgumentParser()
    parser.add_argument("config_path", type=str, default="config.toml")
    # grab a config filepath.
    args = parser.parse_args()
    config = MesospimConfig(args.config_path)
    # Generate a plot for the first laser.
    first_laser_wavelength = config.laser_wavelengths[0]
    t, voltages_t = generate_waveforms(config, first_laser_wavelength)
    # plot the waveforms.
    plot_waveforms_to_pdf(config, t, voltages_t, first_laser_wavelength,
                          f"{first_laser_wavelength}nm_active_plot.pdf")
