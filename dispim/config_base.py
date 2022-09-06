"""TOML wrapper that enables edits, reloads, and manages derived params."""

import tomlkit
import logging
from pathlib import Path
from typing import Union


class TomlConfig:

    def __init__(self, toml_filepath: Union[str, None] = None,
                 config_template: Union[dict, None] = None):
        """Init.

            :param toml_filepath: Optional location of the config if we are
                loading one from file.
            :param config_template: Optional dict with the same key structure
                as the TOML file.
            If specified, (1) a loaded TOML file will be compared with the
            keys in this template and (2) a TOML file can be generated from the
            template.
        """
        self.cfg = None
        self.path = Path(toml_filepath) if toml_filepath else None
        self.template = config_template
        self.log = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        if self.path and self.path.exists():
            self.log.info(f"Loading: {self.path.resolve()}")
            self.load(self.path)
        else:
            self.log.warning("No configuration was specified to load.")
        self.doc_name = self.path.name

    def create_from_template(self):
        """Create a config from a template if one was specified on __init__."""
        if self.template is None:
            raise ValueError("Error: No template was specified from which to "
                             "create the configuration.")
        # This will destroy anything that we loaded.
        self.cfg = tomlkit.parse(tomlkit.dumps(self.template))

    def load(self, filepath: Path = None):
        """Load a config from file specified in filepath or __init__."""
        if filepath:
            self.path = filepath
        assert (self.path.is_file() and self.path.exists()), \
            f"Error: config does not exist at provided filepath: {self.path}."
        with open(filepath, 'r') as toml_file:
            self.cfg = tomlkit.load(toml_file)
        self.doc_name = self.path.name
        if not self.template:
            return
        # TODO: template comparison. possibly with deepdiff.
        #  https://github.com/seperman/deepdiff

    def reload(self):
        """Reload the config from the file we loaded the config from.

        Take all the new changes.
        """
        # This will error out if the config never existed in the first place.
        self.load(self.path)

    def save(self, filepath: str = None, overwrite: bool = True):
        """Save config to specified file, or overwrite if no file specified.

        :param filepath: can be a path to a folder or a file.
            If folder, we use the original filename (or default filename
            if filename never specified).
            If file, we use the specified filename.
            If no path specified, we overwrite unless flagged not to do so.
        :param overwrite: bool to indicate if we overwrite an existing file.
            Defaults to True so that we can save() over a previous file.
        """
        # if filepath unspecified, overwrite the original.
        write_path = Path(filepath) if filepath else self.path
        # if file name is unspecified, use the original.
        if write_path.is_dir():
            write_path = write_path / self.doc_name
        with write_path.open("w") as f:
            f.write(tomlkit.dumps(self.cfg))
