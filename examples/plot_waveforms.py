#!/usr/bin/env python3
"""Waveform Generator to output to pdf."""

from dispim.dispim_config import DispimConfig
from dispim.compute_waveforms import generate_waveforms, plot_waveforms_to_pdf
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--config_path", type=str, default="./config.toml")

    args = parser.parse_args()
    config = DispimConfig(toml_filepath=args.config_path)
    for wavelen in config.laser_wavelengths:
        print(f"Plotting {wavelen}[nm] waveforms")
        t, voltages_t = generate_waveforms(config, wavelen)
        plot_waveforms_to_pdf(config, t, voltages_t, wavelen,
                              f"{wavelen}nm_active_wavelen_plot.pdf")
