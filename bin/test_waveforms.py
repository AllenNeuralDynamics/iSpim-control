from argparse import ArgumentParser
from dispim.dispim_config import DispimConfig
from dispim.compute_waveforms import generate_waveforms, plot_waveforms_to_pdf
from dispim.devices.ni import WaveformHardware
import time

if __name__ == "__main__":
    # Argparse for a config file.
    parser = ArgumentParser()
    parser.add_argument("--config_path", type=str, default="config.toml")
    parser.add_argument("--active_wavelength", type=int, default=488)
    # grab a config filepath.
    args = parser.parse_args()
    config = DispimConfig(args.config_path)
    # Generate a plot for the active laser.
    t_, voltages_of_t_ = generate_waveforms(config, args.active_wavelength)
    # plot the waveforms.
    plot_waveforms_to_pdf(config, t_, voltages_of_t_, args.active_wavelength,
                          f"{args.active_wavelength}nm_active_plot.pdf")

    # Prepare NIDAQ
    ni = WaveformHardware(config.daq_obj_kwds['dev_name'], config.daq_obj_kwds['dev_name'], config.daq_obj_kwds['update_frequency_hz'])
    ni.configure(config.get_daq_cycle_time(), config.daq_ao_names_to_channels)
    ni.assign_waveforms(voltages_of_t_)
    ni.start()
    time.sleep(10)
    ni.stop()
    ni.close()