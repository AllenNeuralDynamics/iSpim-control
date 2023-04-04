#!/usr/bin/env python3
"""Generate waveforms from config params. If standalone, save a graph only."""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from scipy.signal import sawtooth
from ispim.ispim_config import IspimConfig


# TODO: cfg should be able to lookup config params sensibly (like with a string for a key)
def generate_waveforms(cfg: IspimConfig, active_wavelengths: list):
    """return a np nd array with the correct waveforms.
    The DAQ outputs all waveforms to control: the etls, the galvos,
    the cameras, and the laser(s).
    The NI card is triggered by encoder pulses from the stage.
    :param cfg: Mesospim configuration object.
    :param active_wavelength: laser wavelength that will be turned on.
    """
    # initialize empty output voltage array for all interleaved channels
    voltages_out = np.array([]).reshape(cfg.daq_used_channels, 0)
    daq_cycle_time = cfg.get_daq_cycle_time()
    delay_samples = cfg.get_delay_samples()
    exposure_samples = cfg.get_exposure_samples()
    period_samples = cfg.get_period_samples()
    daq_cycle_samples = cfg.get_daq_cycle_samples()
    active_wavelengths.sort()

    for ch in active_wavelengths:
        # Create wavelength-dependent constants
        active_laser_specs = cfg.laser_specs[str(ch)]
        etl_offset = active_laser_specs['etl']['offset']
        etl_amplitude = active_laser_specs['etl']['amplitude']
        galvo_x_offset = active_laser_specs['galvo']['x_offset']
        galvo_y_offset = active_laser_specs['galvo']['y_offset']
        galvo_x_amplitude = active_laser_specs['galvo']['x_amplitude']
        galvo_y_amplitude = active_laser_specs['galvo']['y_amplitude']


        time_samples = np.linspace(0, 2*pi, period_samples)

        # Get peaks of previous sawtooth to connect the waveforms
        previous_galvos = [0,0] if active_wavelengths.index(ch) ==0 \
            else voltages_out[:2,((daq_cycle_samples*active_wavelengths.index(ch))-(daq_cycle_samples-exposure_samples))-1]

        galvo_y, galvo_x = \
                galvo_waveforms(galvo_x_amplitude, galvo_x_offset,
                            galvo_y_amplitude, galvo_y_offset,
                            delay_samples, time_samples, exposure_samples,
                            period_samples, daq_cycle_samples, previous_galvos)

        previous_etl = [0] if active_wavelengths.index(ch) ==0 \
            else voltages_out[2,((daq_cycle_samples*active_wavelengths.index(ch))-(daq_cycle_samples-exposure_samples))-1]

        etl = etl_waveforms(etl_amplitude, etl_offset,
                                  exposure_samples, daq_cycle_samples, previous_etl)
        camera_left = \
            camera_waveforms(exposure_samples, daq_cycle_samples)
        # laser signals arrive in dict, keyed by wavelength in string form.
        laser_signals_dict =\
            laser_waveforms(cfg.laser_specs, ch,
                            exposure_samples, daq_cycle_samples)
        # organize signals by name.
        waveforms = \
        {
            'galvo_y': galvo_y,
            'galvo_x': galvo_x,
            'etl': etl,
            'camera': camera_left,
        }
        waveforms.update(laser_signals_dict)

        # initialize empty output voltage array for single channel
        voltages_t = np.zeros((cfg.daq_used_channels, etl.size))

        # Populate all waveforms in the order the NI card will create them.
        for index, (name, _) in enumerate(cfg.daq_ao_names_to_channels.items()):
            voltages_t[index] = waveforms[name]

        # concatenate and add to the growing output voltage matrix along samples axis

        if active_wavelengths.index(ch) != 0:
            i = active_wavelengths.index(ch)
            voltages_out[:,(daq_cycle_samples*i)-(daq_cycle_samples-exposure_samples):daq_cycle_samples*i] = \
                voltages_t[:,0:daq_cycle_samples-exposure_samples]
            voltages_out = np.concatenate((voltages_out, voltages_t[:,daq_cycle_samples-exposure_samples:]), axis = 1)


        else:

            voltages_out = np.concatenate((voltages_out, voltages_t[:,daq_cycle_samples-exposure_samples:]), axis=1)


    t = np.linspace(0, len(active_wavelengths) * daq_cycle_time, len(active_wavelengths) * daq_cycle_samples,
                    endpoint=False)
    return t, voltages_out


def galvo_waveforms(galvo_x_amplitude, galvo_x_offset,
                    galvo_y_amplitude, galvo_y_offset,
                    delay_samples, time_samples, exposure_samples,
                    period_samples, daq_cycle_samples, previous_peak):

    """Generate galvo waveforms."""
    # Generate relevant galvo signal time chunks
    # Sawtooth duty cycle is adjusted to snapback slowly over the specified rest time
    # x-axis galvos correct for MEMs mirror bow artifact. are quadratic with the y-axis galvo with some scaling amplitude and offset.
    galvo_y = galvo_y_amplitude*sawtooth(time_samples, width=exposure_samples/period_samples) + galvo_y_offset

    # galvo x signal depends on its y signal value.
    galvo_x = abs((galvo_y - galvo_y_offset)**2)*galvo_x_amplitude + galvo_x_offset
    galvo_x[exposure_samples:daq_cycle_samples] = galvo_x[0]  # constant value

    # Adding linearly snapback to beginning of waveform for previous waveform
    start_pos = previous_peak[0]
    end_pos = galvo_y_offset - galvo_y_amplitude
    snap_back = np.linspace(start_pos, end_pos, daq_cycle_samples - exposure_samples)
    galvo_y = np.concatenate((snap_back, galvo_y))

    start_pos = previous_peak[1]
    end_pos = galvo_x[0]
    snap_back = np.linspace(start_pos, end_pos, daq_cycle_samples-exposure_samples)
    galvo_x = np.concatenate((snap_back,galvo_x))

    return galvo_y, galvo_x

def etl_waveforms(etl_amplitude, etl_offset,
                  exposure_samples, daq_cycle_samples, previous_etl):
    """Generate etl waveforms."""
    # ETLs are not actually and are held at DC voltage to correct for axial chromatic shifts.
    etl = etl_amplitude * np.ones(daq_cycle_samples) \
                + etl_offset

    # snap to next etl position as soon as possible
    snapback = np.ones(daq_cycle_samples-exposure_samples)*previous_etl
    etl = np.concatenate((snapback, etl))

    return etl


def camera_waveforms(exposure_samples, daq_cycle_samples):
    """Generate camera waveforms."""
    # Cameras are triggered with some offset specified as % of total exposure time. TODO make this a specified time... not %.
    # Each camera is additionally delay by some time specified by delay time.

    camera_left = np.zeros(daq_cycle_samples)
    camera_left[0:exposure_samples] = 5.0
    snapback = np.zeros(daq_cycle_samples - exposure_samples) # Fake snapback to match other waveforms
    camera_left = np.concatenate((snapback,camera_left))

    return camera_left


def laser_waveforms(laser_specs, active_wavelen: int,
                    exposure_samples, daq_cycle_samples):
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
        laser_t = disable_voltage*np.ones(daq_cycle_samples)
        laser_t[0:exposure_samples] = enable_voltage

        snapback = np.zeros(daq_cycle_samples - exposure_samples)  # Fake snapback to match other waveforms
        laser_t = np.concatenate((snapback, laser_t))

        lasers_t[ao_name] = laser_t
    return lasers_t


def plot_waveforms_to_pdf(cfg: IspimConfig, t: np.array,
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
    from ispim_config import IsipimConfig
    # Argparse for a config file.
    parser = ArgumentParser()
    parser.add_argument("config_path", type=str, default="config.toml")
    parser.add_argument("active_wavelength", type=int, default=488)
    # grab a config filepath.
    args = parser.parse_args()
    config = IspimConfig(args.config_path)
    # Generate a plot for the active laser.
    t_, voltages_of_t_ = generate_waveforms(config, args.active_wavelength)
    # plot the waveforms.
    plot_waveforms_to_pdf(config, t_, voltages_of_t_, args.active_wavelength,
                          f"{args.active_wavelength}nm_active_plot.pdf")