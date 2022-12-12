#!/usr/bin/env python3
"""Command line utility for generating a default dispim configuration."""

from dispim.dispim_config import DispimConfig
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Create a default dispim configuration from the template.")
    parser.add_argument("--output_path", "-o", type=str,
                        default="./config.toml")

    args = parser.parse_args()
    config = DispimConfig(args.output_path, create=True)
    config.save()

