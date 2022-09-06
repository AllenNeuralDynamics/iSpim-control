"""Spim Base class."""

import pkg_resources
import logging
import shutil
import sys
from coloredlogs import ColoredFormatter
from contextlib import contextmanager
from datetime import date
from git import Repo
from pathlib import Path


class Spim:

    def __init__(self, config_filepath: str,
                 log_filename: str = 'debug.log',
                 console_output: bool = True,
                 color_console_output: bool = False,
                 console_output_level: str = 'info',
                 simulated: bool = False):
        """Read config file. Create Mesopim components according to config."""
        # If simulated, no physical connections to hardware should be required.
        # Simulation behavior should be pushed as far down the object
        # hierarchy as possible.
        self.simulated = simulated
        # Setup logging.
        # Save console output to print/not-print imaging progress.
        self.console_output = console_output
        # We want the name of the package here since logger hierarchy
        # depends on module structure.
        self.log = logging.getLogger(__package__)
        # logger level must be set to the lowest level of any handler.
        self.log.setLevel(logging.DEBUG)
        # Create log handlers to dispatch:
        # - DEBUG level and above to write to a file called debug.log.
        # - User-specified level and above to print to console if specified.
        fmt = '%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s'
        fmt = "[SIM] " + fmt if self.simulated else fmt
        datefmt = '%Y-%m-%d,%H:%M:%S'
        self.log_format = logging.Formatter(fmt=fmt, datefmt=datefmt)
        self.log_handlers = []
        debug_filepath = Path(log_filename)
        self.log_handlers.append(logging.FileHandler(debug_filepath, 'w'))
        self.log_handlers[-1].setLevel(logging.DEBUG)
        self.log_handlers[-1].setFormatter(self.log_format)
        if self.console_output:
            self.log_handlers.append(logging.StreamHandler(sys.stdout))
            self.log_handlers[-1].setLevel(console_output_level)
            if color_console_output:
                colored_formatter = ColoredFormatter(fmt=fmt, datefmt=datefmt)
                self.log_handlers[-1].setFormatter(colored_formatter)
            else:
                self.log_handlers[-1].setFormatter(self.log_format)
        for handler in self.log_handlers:
            self.log.addHandler(handler)

        # Location where images will be saved to from a config-based run.
        # We don't know this in advance since we create the folder at runtime.
        self.img_storage_dir = None

        # Config. This will get defined in the child class.
        self.cfg = None

    def log_git_hashes(self):
        """Log the git hashes of this project and all packages."""
        # Iterate through this pkg's required packages and log all git hashes.
        # Warn if they have been changed.
        env = {str(ws): ws.module_path for ws in pkg_resources.working_set if
               ws.module_path and  # This can be None.
               not ws.module_path.endswith(('site-packages', 'dist-packages'))}
        for pkg_name, env_path in env.items():
            repo = Repo(env_path)
            self.log.debug(f"{pkg_name} on branch {repo.active_branch} at "
                           f"{repo.head.object.hexsha}")
            if repo.is_dirty():
                self.log.error(f"{pkg_name} has uncommitted changes.")

    def log_runtime_estimate(self):
        """Log how much time the configured imaging run will take."""
        raise NotImplementedError

    @contextmanager
    def log_to_file(self, log_filepath: str):
        """Log to a file for the duration of a function's execution."""
        try:
            log_handler = logging.FileHandler(log_filepath, 'w')
            log_handler.setLevel(logging.DEBUG)
            log_handler.setFormatter(self.log_format)
            self.log.addHandler(log_handler)
            yield
        finally:
            log_handler.close()
            self.log.removeHandler(log_handler)

    def run(self, overwrite: bool = False):
        """Collect data according to config; populate dest folder with data.

        :param overwrite: bool indicating if we want to overwrite any existing
            data if the output folder already exists. False by default.
        """
        # Create output folder and folder for storing images.
        output_folder = \
            self.cfg.ext_storage_dir / Path(self.cfg.subject_id + "-ID_" +
                                            date.today().strftime("%Y_%m_%d"))
        print(f"external storage dir is: {self.cfg.ext_storage_dir.resolve()}")
        if output_folder.exists() and not overwrite:
            self.log.error(f"Output folder {output_folder.absolute()} exists, "
                           "This function must be rerun with overwrite=True.")
            raise
        self.img_storage_dir = output_folder / Path("micr/")
        self.log.info(f"Creating datset folder in: {output_folder.absolute()}")
        self.img_storage_dir.mkdir(parents=True, exist_ok=overwrite)
        # Save the config file we will run.
        self.cfg.save(output_folder, overwrite=overwrite)
        # Log to a file for the duration of this function's execution.
        imaging_log_filepath = Path("imaging_log.log")  # name should be a constant.
        with self.log_to_file(imaging_log_filepath):
            self.log_git_hashes()
            self.run_from_config()
        # Copy the log file we just made to the folder with our data.
        # Note: shutil can't overwrite, so we must delete any prior imaging log
        #   in the destination folder if we are overwriting.
        imaging_log_dest = output_folder/Path(imaging_log_filepath.name)
        if overwrite and imaging_log_dest.exists():
            imaging_log_dest.unlink()
        # We must use shutil because we may be moving files across disks.
        shutil.move(str(imaging_log_filepath), str(output_folder))

    def run_from_config(self):
        raise NotImplementedError("Child class must implement this function.")

    def livestream(self):
        raise NotImplementedError

    def reload_config(self):
        """Reload the toml file."""
        self.cfg.reload()

    def close(self):
        """Safely close all open hardware connections."""
        # Most of the action here should be implemented in a child class that
        # calls super().close() at the very end.
        self.log.info("Ending log.")
        for handler in self.log_handlers:
            handler.close()
