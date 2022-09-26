#!/usr/bin/env python3
"""Generate waveforms from config params. If standalone, save a graph only."""
import np as np
from np import pi
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
    etl_left_offset = active_laser_specs['etl_left']['offset']
    etl_right_offset = active_laser_specs['etl_right']['offset']
    etl_left_amplitude = active_laser_specs['etl_left']['amplitude']
    etl_right_amplitude = active_laser_specs['etl_right']['amplitude']
    galvo_x_left_offset = active_laser_specs['galvo_x_left']['offset']
    galvo_x_right_offset = active_laser_specs['galvo_x_right']['offset']
    galvo_y_left_offset = active_laser_specs['galvo_y_left']['offset']
    galvo_y_right_offset = active_laser_specs['galvo_y_right']['offset']
    galvo_x_left_amplitude = active_laser_specs['galvo_x_left']['amplitude']
    galvo_x_right_amplitude = active_laser_specs['galvo_x_right']['amplitude']
    galvo_y_left_amplitude = active_laser_specs['galvo_y_left']['amplitude']
    galvo_y_right_amplitude = active_laser_specs['galvo_y_right']['amplitude']
    camera_left_offset = cfg.camera_left_offset
    camera_right_offset = cfg.camera_right_offset

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
    laser = laser_waveforms(camera_left_offset, camera_right_offset,
                            delay_samples, exposure_samples, daq_cycle_samples)
    # Create mapping to stuff waveforms in the same order that we create them
    # on the NI card.
    waveforms = \
    {
        'galvo_y_left': galvo_y_left,
        'galvo_y_right': galvo_y_right,
        'galvo_x_left': galvo_x_left,
        'galvo_x_right': galvo_x_right,
        'etl_left': etl_left,
        'etl_right': etl_right,
        'camera_left': camera_left,
        'camera_right': camera_right
    }
    for index, name, _ in enumerate(cfg.daq_ao_names_to_channels.items()):
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
    galvo_x_left_t[delay_samples:delay_samples+period_samples] = galvo_x_left
    galvo_x_left_t[period_samples::] = galvo_x_left[-1]  # constant value
    galvo_x_right_t[0:period_samples] = galvo_x_right
    galvo_x_right_t[period_samples::] = galvo_x_right[0]  # constant value
    galvo_y_left_t[delay_samples:delay_samples+period_samples] = galvo_y_left
    galvo_y_left[period_samples::] = galvo_y_left[-1]  # constant value
    galvo_y_right[0:period_samples] = galvo_y_right
    galvo_y_right[period_samples::] = galvo_y_right[0]  # constant value

    return galvo_y_left_t, galvo_y_right_t, galvo_x_left_t, galvo_x_right_t


def etl_waveforms(etl_left_amplitude, etl_left_offset,
                  etl_right_amplitude, etl_right_offset,
                  daq_cycle_samples):
    """Generate etl waveforms."""
    # ETLs are not actually and are held at DC voltage to correct for axial chromatic shifts.
    etl_left = etl_left_amplitude * np.ones(len(daq_cycle_samples)) \
               + etl_left_offset
    etl_right = etl_right_amplitude * np.ones(len(daq_cycle_samples)) \
                + etl_right_offset

    return etl_left, etl_right


def camera_waveforms(camera_left_offset, camera_right_offset,
                     delay_samples, exposure_samples, daq_cycle_samples):
    """Generate camera waveforms."""
    # Cameras are triggered with some offset specified as % of total exposure time. TODO make this a specified time... not %.
    # Each camera is additionally delay by some time specified by delay time.
    camera_right = np.zeros((1,daq_cycle_samples))
    camera_right[int(exposure_samples*camera_right_offset):int(exposure_samples*camera_right_offset) + exposure_samples] = 5.0
    camera_left = np.zeros((1,daq_cycle_samples))
    camera_left[delay_samples + int(exposure_samples*camera_left_offset):delay_samples + int(exposure_samples*camera_left_offset) + exposure_samples] = 5.0

    return camera_left, camera_right


def laser_waveforms(camera_left_offset, camera_right_offset,
                    delay_samples, exposure_samples, daq_cycle_samples):
    """Generate laser waveforms."""
    # Lasers are triggered and strobed only when the camera is exposing.
    laser = np.zeros((1,daq_cycle_samples))
    laser[int(exposure_samples*camera_right_offset):delay_samples + int(exposure_samples*camera_left_offset) + exposure_samples] = 5.0

    return laser


def plot_waveforms_to_pdf(cfg: DispimConfig, t: np.array,
                          voltages_t: np.array, active_wavelength: int,
                          filename: str = "plot.pdf"):

        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 7))

        axes.set_title("One Image Capture Sequence.")
        axes.plot(t, voltages_t[cfg.daq_ao_names_to_channels['galvo_x_left']], label="galvo_x_left")
        axes.plot(t, voltages_t[cfg.daq_ao_names_to_channels['galvo_x_right']], label="galvo_x_right")
        axes.plot(t, voltages_t[cfg.daq_ao_names_to_channels['galvo_y_left']], label="galvo_y_left")
        axes.plot(t, voltages_t[cfg.daq_ao_names_to_channels['galvo_y_right']], label="galvo_y_right")
        axes.plot(t, voltages_t[cfg.daq_ao_names_to_channels['etl_left']], label="etl_left")
        axes.plot(t, voltages_t[cfg.daq_ao_names_to_channels['etl_right']], label="etl_right")
        axes.plot(t, voltages_t[cfg.daq_ao_names_to_channels['camera_left']], label="camera_left")
        axes.plot(t, voltages_t[cfg.daq_ao_names_to_channels['camera_right']], label="camera_right")
        axes.plot(t, voltages_t[cfg.daq_ao_names_to_channels['laser']], label=str(active_wavelength))
        axes.set_xlabel("time [s]")
        axes.set_ylabel("amplitude [V]")
        axes.set_ylim(0,5)
        axes.legend(loc="center")

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
    # grab a config filepath.
    args = parser.parse_args()
    config = DispimConfig(args.config_path)
    # Generate a plot for the active laser. FIX HERE TO GRAB ACTIVE WAVELENGTH
    active_wavelength = config.laser_wavelengths[0] # FIX HERE TO GRAB ACTIVE WAVELENGTH
    t, voltages_t = generate_waveforms(config, active_wavelength)
    # plot the waveforms.
    plot_waveforms_to_pdf(config, t, voltages_t, active_wavelength,
                          f"{active_wavelength}nm_active_plot.pdf")