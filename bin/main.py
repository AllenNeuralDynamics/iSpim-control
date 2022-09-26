#!/usr/bin/env python3
"""example main script to launch the mesospim."""

from dispim.dispim import Dispim
from coloredlogs import ColoredFormatter
import ctypes
import logging
import argparse
import os
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, default=None)
    parser.add_argument("--log_level", type=str, default="INFO",
                        choices=["INFO", "DEBUG"])
    parser.add_argument("--simulated", default=False, action="store_true",
                        help="Simulate hardware device connections.")
    parser.add_argument("--console_output", default=True,
                        help="whether or not to print to the console.")
    # Note: colored console output is buggy on Windows.
    parser.add_argument("--color_console_output", action="store_true",
                        default=False if os.name == 'nt' else True)

    args = parser.parse_args()
    # Check if we didn't supply a config file and populate a safe guess.
    if not args.config:
        if args.simulated:
            args.config = "./sim_config.toml"
        else:
            args.config = "./config.toml"

    # Setup logging.
    # Create log handlers to dispatch:
    # - DEBUG level and above to write to a file called debug.log.
    # - User-specified level and above to print to console if specified.
    logger = logging.getLogger()  # get the root logger.
    # logger level must be set to the lowest level of any handler.
    logger.setLevel(logging.DEBUG)
    fmt = '%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s'
    fmt = "[SIM] " + fmt if args.simulated else fmt
    datefmt = '%Y-%m-%d,%H:%M:%S'
    log_format = logging.Formatter(fmt=fmt, datefmt=datefmt)
    log_handlers = []
    #debug_filepath = Path(log_filename)
    #self.log_handlers.append(logging.FileHandler(debug_filepath, 'w'))
    #log_handlers[-1].setLevel(logging.DEBUG)
    #log_handlers[-1].setFormatter(log_format)
    if args.console_output:
        log_handlers.append(logging.StreamHandler(sys.stdout))
        log_handlers[-1].setLevel(args.log_level)
        if args.color_console_output:
            colored_formatter = ColoredFormatter(fmt=fmt, datefmt=datefmt)
            log_handlers[-1].setFormatter(colored_formatter)
        else:
            log_handlers[-1].setFormatter(log_format)
    for handler in log_handlers:
        logger.addHandler(handler)

    # Windows-based console needs to accept colored logs if running with color.
    if os.name == 'nt' and args.color_console_output:
        kernel32 = ctypes.windll.kernel32
        kernel32.SetConsoleMode(kernel32.GetStdHandle(-11), 7)

    instrument = Dispim(config_filepath=args.config,
                          #console_output_level=args.log_level,
                          #color_console_output=args.color_console_output,
                          simulated=args.simulated)
    try:
        #from inpromptu import Inpromptu
        #Inpromptu(instrument).cmdloop()
        instrument.run(overwrite=args.simulated)
    except KeyboardInterrupt:
        pass
    finally:
        instrument.close()

if __name__ == '__main__':
    main()
