#!/usr/bin/env python3
"""Generate waveforms from config params. If standalone, save a graph only."""
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from scipy.signal import square, sawtooth
from .mesospim_config import MesospimConfig

# TODO: cfg should be able to lookup config params sensibly (like with a string for a key)
def generate_waveforms(cfg: MesospimConfig, active_wavelength: int):
    """return a numpy nd array with the correct waveforms.
    The DAQ outputs all waveforms to control: the etls, the galvos,
    the cameras, and the laser(s).
    The NI card is triggered by encoder pulses from the stage.
    :param cfg: Mesospim configuration object.
    :param active_wavelength: laser wavelength that will be turned on.
    """
    # Create wavelength-dependent constants
    self.active_laser_specs = cfg.laser_specs[str(active_wavelength)]
    self.etl_left_offset = active_laser_specs['etl_left']['offset']
    self.etl_right_offset = active_laser_specs['etl_right']['offset']
    self.etl_left_amplitude = active_laser_specs['etl_left']['amplitude']
    self.etl_right_amplitude = active_laser_specs['etl_right']['amplitude']
    self.galvo_x_left_offset = active_laser_specs['galvo_x_left']['offset']
    self.galvo_x_right_offset = active_laser_specs['galvo_x_right']['offset']
    self.galvo_y_left_offset = active_laser_specs['galvo_y_left']['offset']
    self.galvo_y_right_offset = active_laser_specs['galvo_y_right']['offset']
    self.galvo_x_left_amplitude = active_laser_specs['galvo_x_left']['amplitude']
    self.galvo_x_right_amplitude = active_laser_specs['galvo_x_right']['amplitude']
    self.galvo_y_left_amplitude = active_laser_specs['galvo_y_left']['amplitude']
    self.galvo_y_right_amplitude = active_laser_specs['galvo_y_right']['amplitude']
    self.camera_left_offset = cfg.camera_left_offset
    self.camera_right_offset = cfg.camera_right_offset

    self.delay_time = cfg.delay_time
    self.rest_time = cfg.rest_time
    self.exposure_time = cfg.exposure_time
    self.period_time = cfg.get_period_time()
    self.daq_cycle_time = cfg.get_daq_cycle_time()

    self.delay_samples = cfg.get_delay_samples()
    self.rest_samples = cfg.get_rest_samples()
    self.exposure_samples = cfg.get_exposure_samples()
    self.period_samples = cfg.get_period_samples()
    self.daq_cycle_samples = cfg.get_daq_cycle_samples()

    self.time_samples = numpy.linspace(0, 2*math.pi, period_samples)

    # Create table that holds an entire cycle's worth of pts.
    # External trigger signal is the last value. Aside from that,
    # the signal order in voltages_t doesn't matter provide that it is
    # consistent with the signal order that the WaveformGenerator creates
    # tasks. We can ensure this by iterating through the same data structure
    # in the cfg as we do when we create tasks with the WaveformGenerator.

    voltages_t = np.zeros((cfg.daq_used_channels, daq_cycle_samples))
    t = np.linspace(0, daq_cycle_time, daq_num_samples, endpoint=False)

    galvo_y_left, galvo_y_right, galvo_x_left, galvo_x_right = galvo_waveforms()
    etl_left, etl_right = etl_waveforms()
    camera_left, camera_right = camera_waveforms()
    laser = laser_waveforms()

    # Add left y galvo waveforms
    voltages_t[cfg.daq_ao_names_to_channels['galvo_y_left'], self.delay_samples:self.delay_samples+self.period_samples] = galvo_y_left
    voltages_t[cfg.daq_ao_names_to_channels['galvo_y_left'], self.period_samples::] = galvo_y_left[-1]
    # Add right y galvo waveforms
    voltages_t[cfg.daq_ao_names_to_channels['galvo_y_right'], 0:self.period_samples] =  galvo_y_right
    voltages_t[cfg.daq_ao_names_to_channels['galvo_y_right'], self.period_samples::] = galvo_y_right[0]
    # Add left x galvo waveforms
    voltages_t[cfg.daq_ao_names_to_channels['galvo_x_left'], self.delay_samples:self.delay_samples+self.period_samples] = galvo_x_left
    voltages_t[cfg.daq_ao_names_to_channels['galvo_x_left'], self.period_samples::] = galvo_x_left[-1]
    # Add right x galvo waveforms
    voltages_t[cfg.daq_ao_names_to_channels['galvo_x_right'], 0:self.period_samples] =  galvo_x_right
    voltages_t[cfg.daq_ao_names_to_channels['galvo_x_right'], self.period_samples::] = galvo_x_right[0]
    # Add etl  waveforms
    voltages_t[cfg.daq_ao_names_to_channels['etl_left'], :] = etl_left
    voltages_t[cfg.daq_ao_names_to_channels['etl_right'], :] = etl_right
    # Add camera waveforms
    voltages_t[cfg.daq_ao_names_to_channels['camera_left'], :] = camera_left
    voltages_t[cfg.daq_ao_names_to_channels['camera_right'], :] = camera_right
    # Add active laser waveforms
    voltages_t[cfg.daq_ao_names_to_channels[str(active_wavelength)], :] = laser

    return t, voltages_t

def galvo_waveforms():
    """Generate galvo waveforms."""
    # Sawtooth duty cycle is adjusted to snapback slowly over the specified rest time
    # x-axis galvos correct for MEMs mirror bow artifact. are quadratic with the y-axis galvo with some scaling amplitude and offset.
    galvo_y_left = -self.galvo_y_left_amplitude*signal.sawtooth(self.time_samples, width=self.exposure_samples*1.0/self.period_samples) + self.galvo_y_left_offset
    galvo_y_right = self.galvo_y_right_amplitude*signal.sawtooth(self.time_samples, width=self.exposure_samples*1.0/self.period_samples) + self.galvo_y_right_offset
    galvo_x_left = abs((self.galvo_y_left_amplitude - self.galvo_y_left_offset)**2)*self.galvo_x_left_amplitude + self.galvo_x_left_offset
    galvo_x_right = abs((self.galvo_y_right_amplitude - self.galvo_y_right_offset)**2)*self.galvo_x_right_amplitude + self.galvo_x_right_offset

    return galvo_y_left, galvo_y_right, galvo_x_left, galvo_x_right

def etl_waveforms():
    """Generate etl waveforms."""
    # ETLs are not actually and are held at DC voltage to correct for axial chromatic shifts.
    etl_left = self.etl_amplitude_left
    etl_right = self.etl_amplitude_right

    return etl_left, etl_right

def camera_waveforms():
    """Generate camera waveforms."""
    # Cameras are triggered with some offset specified as % of total exposure time. TODO make this a specified time... not %.
    # Each camera is additionally delay by some time specified by delay time.
    camera_right = numpy.zeros((1,self.daq_cycle_samples))
    camera_right[int(self.exposure_samples*self.camera_right_offset):int(self.exposure_samples*self.camera_right_offset) + self.exposure_samples] = 5.0
    camera_left = numpy.zeros((1,self.daq_cycle_samples))
    camera_left[self.delay_samples + int(self.exposure_samples*self.camera_left_offset):self.delay_samples + int(self.exposure_samples*self.camera_left_offset) + self.exposure_samples] = 5.0

    return camera_left, camera_right

def laser_waveforms():
    """Generate laser waveforms."""
    # Lasers are triggered and strobed only when the camera is exposing.
    laser = numpy.zeros((1,self.daq_cycle_samples))
    laser[int(self.exposure_samples*self.camera_right_offset):self.delay_samples + int(self.exposure_samples*self.camera_left_offset) + self.exposure_samples] = 5.0

    return laser

def plot_waveforms_to_pdf(cfg: MesospimConfig, t: np.array,
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
    from .config_base import MesospimConfig
    # Argparse for a config file.
    parser = ArgumentParser()
    parser.add_argument("config_path", type=str, default="config.toml")
    # grab a config filepath.
    args = parser.parse_args()
    config = MesospimConfig(args.config_path)
    # Generate a plot for the active laser. FIX HERE TO GRAB ACTIVE WAVELENGTH
    active_wavelength = config.laser_wavelengths[0] # FIX HERE TO GRAB ACTIVE WAVELENGTH
    t, voltages_t = generate_waveforms(config, active_wavelength)
    # plot the waveforms.
    plot_waveforms_to_pdf(config, t, voltages_t, active_wavelength,
                          f"{active_wavelength}nm_active_plot.pdf")