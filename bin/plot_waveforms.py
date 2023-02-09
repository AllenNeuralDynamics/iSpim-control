from argparse import ArgumentParser
from dispim.dispim_config import DispimConfig
from dispim.compute_waveforms import generate_waveforms, plot_waveforms_to_pdf

if __name__ == "__main__":
    # Argparse for a config file.
    parser = ArgumentParser()
    parser.add_argument("--config_path", type=str, default="config.toml")
    parser.add_argument("--active_wavelength", type=int, default=[488, 561])
    # grab a config filepath.
    args = parser.parse_args()
    config = DispimConfig(args.config_path)
    # Generate a plot for the active laser.
    t_, voltages_of_t_ = generate_waveforms(config, args.active_wavelength)
    # plot the waveforms.
    plot_waveforms_to_pdf(config, t_, voltages_of_t_, args.active_wavelength,
                          f"{args.active_wavelength}nm_active_plot.pdf")