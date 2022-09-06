#!/usr/bin/env python3
"""Command line utility for generating a blank mesospim configuration."""
import time

from mesospim.mesospim_config import MesospimConfig
import argparse
#import argcomplete

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_path", type=str, default="./config.toml")
    #argcomplete.autocomplete(parser)

    args = parser.parse_args()
    config = MesospimConfig(args.output_path)
    #config.create_from_template()
    #config.save(args.output_path)

