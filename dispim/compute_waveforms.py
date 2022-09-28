#!/usr/bin/env python3
"""Generate waveforms from config params. If standalone, save a graph only."""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from scipy.signal import sawtooth
from .dispim_config import DispimConfig


# TODO: cfg should be able to lookup config params sensibly (like with a string for a key)
def generate_waveforms(cfg: DispimConfig, active_wavelength: int):
    """return a np nd array with the correct waveforms.
    The DAQ outputs all waveforms to control: the etls, the galvos,
    the cameras, and the laser(s).
    The NI card is triggered by encoder pulses from the stage.
    :param cfg: Mesospim configuration object.
    :param active_wavelength: laser wavelength that will be turned on.
    """
    # Create wavelength-dependent constants
    active_laser_specs = cfg.laser_specs[str(active_wavelength)]
    etl_left_offset = active_laser_specs['etl']['left_offset']
    etl_right_offset = active_laser_specs['etl']['right_offset']
    etl_left_amplitude = active_laser_specs['etl']['left_amplitude']
    etl_right_amplitude = active_laser_specs['etl']['right_amplitude']
    galvo_x_left_offset = active_laser_specs['galvo']['x_left_offset']
    galvo_x_right_offset = active_laser_specs['galvo']['x_right_offset']
    galvo_y_left_offset = active_laser_specs['galvo']['y_left_offset']
    galvo_y_right_offset = active_laser_specs['galvo']['y_right_offset']
    galvo_x_left_amplitude = active_laser_specs['galvo']['x_left_amplitude']
    galvo_x_right_amplitude = active_laser_specs['galvo']['x_right_amplitude']
    galvo_y_left_amplitude = active_laser_specs['galvo']['y_left_amplitude']
    galvo_y_right_amplitude = active_laser_specs['galvo']['y_right_amplitude']
    camera_left_offset = cfg.camera_specs['camera_left']['offset']
    camera_right_offset = cfg.camera_specs['camera_right']['offset']

    daq_cycle_time = cfg.get_daq_cycle_time()
    delay_samples = cfg.get_delay_samples()
    exposure_samples = cfg.get_exposure_samples()
    period_samples = cfg.get_period_samples()
    daq_cycle_samples = cfg.get_daq_cycle_samples()

    time_samples = np.linspace(0, 2*pi, period_samples)

    voltages_t = np.zeros((cfg.daq_used_channels, daq_cycle_samples))
    t = np.linspace(0, daq_cycle_time, daq_cycle_samples, endpoint=False)

    # Create time serieses for each specified signal.
    galvo_y_left, galvo_y_right, galvo_x_left, galvo_x_right = \
        galvo_waveforms(galvo_x_left_amplitude, galvo_x_left_offset,
                        galvo_x_right_amplitude, galvo_x_right_offset,
                        galvo_y_left_amplitude, galvo_y_left_offset,
                        galvo_y_right_amplitude, galvo_y_right_offset,
                        delay_samples, time_samples, exposure_samples,
                        period_samples, daq_cycle_samples)
    etl_left, etl_right = etl_waveforms(etl_left_amplitude, etl_left_offset,
                                        etl_right_amplitude, etl_right_offset,
                                        daq_cycle_samples)
    camera_left, camera_right = \
        camera_waveforms(camera_left_offset, camera_right_offset,
                         delay_samples, exposure_samples, daq_cycle_samples)
    # laser signals arrive in dict, keyed by wavelength in string form.
    laser_signals_dict =\
        laser_waveforms(cfg.laser_specs, active_wavelength,
                        camera_left_offset, camera_right_offset,
                        delay_samples, exposure_samples, daq_cycle_samples)
    # organize signals by name.
    waveforms = \
    {
        'galvo_y_left': galvo_y_left,
        'galvo_y_right': galvo_y_right,
        'galvo_x_left': galvo_x_left,
        'galvo_x_right': galvo_x_right,
        'etl_left': etl_left,
        'etl_right': etl_right,
        'camera_left': camera_left,
        'camera_right': camera_right,
    }
    waveforms.update(laser_signals_dict)

    # Populate all waveforms in the order the NI card will create them.
    for index, (name, _) in enumerate(cfg.daq_ao_names_to_channels.items()):
        voltages_t[index] = waveforms[name]

    return t, voltages_t


def galvo_waveforms(galvo_x_left_amplitude, galvo_x_left_offset,
                    galvo_x_right_amplitude, galvo_x_right_offset,
                    galvo_y_left_amplitude, galvo_y_left_offset,
                    galvo_y_right_amplitude, galvo_y_right_offset,
                    delay_samples, time_samples, exposure_samples,
                    period_samples, daq_cycle_samples):
    """Generate galvo waveforms."""
    # Create full signal as played by the daq.
    galvo_x_left_t = np.zeros(daq_cycle_samples)
    galvo_x_right_t = np.zeros(daq_cycle_samples)
    galvo_y_left_t = np.zeros(daq_cycle_samples)
    galvo_y_right_t = np.zeros(daq_cycle_samples)

    # Generate relevant galvo signal time chunks
    # Sawtooth duty cycle is adjusted to snapback slowly over the specified rest time
    # x-axis galvos correct for MEMs mirror bow artifact. are quadratic with the y-axis galvo with some scaling amplitude and offset.
    galvo_y_left = -galvo_y_left_amplitude*sawtooth(time_samples, width=exposure_samples*1.0/period_samples) + galvo_y_left_offset
    galvo_y_right = galvo_y_right_amplitude*sawtooth(time_samples, width=exposure_samples*1.0/period_samples) + galvo_y_right_offset
    # galvo x signal depends on its y signal value.
    galvo_x_left = abs((galvo_y_left - galvo_y_left_offset)**2)*galvo_x_left_amplitude + galvo_x_left_offset
    galvo_x_right = abs((galvo_y_right - galvo_y_right_offset)**2)*galvo_x_right_amplitude + galvo_x_right_offset

    # stuff galvo signal time chunks in the right location.
    galvo_y_left_t[delay_samples:delay_samples+period_samples] = galvo_y_left
    galvo_y_left_t[0:delay_samples] = galvo_y_left[-1]  # constant value

    galvo_y_right_t[0:period_samples] = galvo_y_right
    galvo_y_right_t[period_samples::] = galvo_y_right[0]  # constant value

    galvo_x_left_t[delay_samples:delay_samples+period_samples] = galvo_x_left
    galvo_x_left_t[0:delay_samples] = galvo_x_left[-1]   # constant value
    galvo_x_left_t[delay_samples+exposure_samples::] = galvo_x_left[-1]  # constant value

    galvo_x_right_t[0:period_samples] = galvo_x_right
    galvo_x_right_t[period_samples::] = galvo_y_right[0]  # constant value
    galvo_x_right_t[exposure_samples::] = galvo_x_right[-1]  # constant value

    return galvo_y_left_t, galvo_y_right_t, galvo_x_left_t, galvo_x_right_t


def etl_waveforms(etl_left_amplitude, etl_left_offset,
                  etl_right_amplitude, etl_right_offset,
                  daq_cycle_samples):
    """Generate etl waveforms."""
    # ETLs are not actually and are held at DC voltage to correct for axial chromatic shifts.
    etl_left = etl_left_amplitude * np.ones(daq_cycle_samples) \
               + etl_left_offset
    etl_right = etl_right_amplitude * np.ones(daq_cycle_samples) \
                + etl_right_offset

    return etl_left, etl_right


def camera_waveforms(camera_left_offset, camera_right_offset,
                     delay_samples, exposure_samples, daq_cycle_samples):
    """Generate camera waveforms."""
    # Cameras are triggered with some offset specified as % of total exposure time. TODO make this a specified time... not %.
    # Each camera is additionally delay by some time specified by delay time.
    camera_right = np.zeros(daq_cycle_samples)
    camera_right[int(exposure_samples*camera_right_offset):int(exposure_samples*camera_right_offset) + exposure_samples] = 5.0
    camera_left = np.zeros(daq_cycle_samples)
    camera_left[delay_samples + int(exposure_samples*camera_left_offset):delay_samples + int(exposure_samples*camera_left_offset) + exposure_samples] = 5.0

    return camera_left, camera_right


def laser_waveforms(laser_specs, active_wavelen: int,
                    camera_left_offset, camera_right_offset,
                    delay_samples, exposure_samples, daq_cycle_samples):
    """Generate multiple laser waveforms with one active signal and the rest
    flatlined.

    :return: dictionary, keyed by string form of laser wavelength in nm.
    """
    lasers_t = {}
    # Generate Laser Signal. Signal should be:
    # a pulse for the active laser but a flatline for inactive lasers.
    for ao_name in [str(nm) for nm in laser_specs]:
        disable_voltage = laser_specs[ao_name]['disable_voltage']
        enable_voltage = laser_specs[ao_name]['disable_voltage']
        # Only enable the active wavelengths.
        if ao_name == str(active_wavelen):
            enable_voltage = laser_specs[ao_name]['enable_voltage']
        # Generate Laser Signal analog time series.
        laser_t = disable_voltage * np.ones((1, daq_cycle_samples))
        laser_t[int(exposure_samples*camera_right_offset):delay_samples + int(exposure_samples*camera_left_offset) + exposure_samples] = enable_voltage
        lasers_t[ao_name] = laser_t
    return lasers_t


def plot_waveforms_to_pdf(cfg: DispimConfig, t: np.array,
                          voltages_t: np.array, active_wavelength: int,
                          filename: str = "plot.pdf"):

        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 7))

        axes.set_title("One Image Capture Sequence.")
        for index, ao_name in enumerate(cfg.daq_ao_names_to_channels.keys()):
            axes.plot(t, voltages_t[index], label=ao_name)
        axes.set_xlabel("time [s]")
        axes.set_ylabel("amplitude [V]")
        # axes.set_ylim(2,3)
        axes.legend(loc="upper right")

        try:
            fig.savefig(filename)
        except OSError as e:
            print("Error: cannot save figure. Another program may be using it.")
            raise e


if __name__ == "__main__":
    from argparse import ArgumentParser
    from .dispim_config import DispimConfig
    # Argparse for a config file.
    parser = ArgumentParser()
    parser.add_argument("config_path", type=str, default="config.toml")
    parser.add_argument("active_wavelength", type=int, default=488)
    # grab a config filepath.
    args = parser.parse_args()
    config = DispimConfig(args.config_path)
    # Generate a plot for the active laser.
    t_, voltages_of_t_ = generate_waveforms(config, args.active_wavelength)
    # plot the waveforms.
    plot_waveforms_to_pdf(config, t_, voltages_of_t_, args.active_wavelength,
                          f"{args.active_wavelength}nm_active_plot.pdf")