"""Mesospim Config object to flatten item access in the TOML-based config."""

from mesospim.config_base import SpimConfig

# A template from which we can generate a blank mesospim config toml file.
TomlTemplate = \
    {
        "imaging_specs":
            {
                "local_storage_directory": ".",
                "external_storage_directory": ".",
                "subject_id": "brain_00",
                "tile_prefix": "tile",
                "tile_overlap_x_percent": 15,
                "tile_overlap_y_percent": 15,
            }
        # TODO: populate remaining fields.
    }


class DispimConfig(SpimConfig):
    """A Dispim Configuration."""

    def __init__(self, toml_filepath: str):
        """Read config file. Warn if not found, but create sensible defaults."""
        super().__init__(toml_filepath, TomlTemplate)

        # Note: these are mutable, so reloading the toml doesn't affect them.
        self.stage_specs = self.cfg['sample_stage_specs']
        self.laser_specs = self.cfg['channel_specs']

        self.tiger_obj_kwds = self.cfg['tiger_controller_driver_kwds']
        self.daq_obj_kwds = self.cfg['daq_driver_kwds']
        # TODO: dispim has 2 filterwheels. We must set the location of both
        #   programmatically.
        #self.filter_wheel_kwds = self.cfg['filter_wheel_kwds']

        # constants.
        # TODO: put any config-related constants here.

    # Getters. These must be functions since values depend on the laser
    # wavelength. Otherwise, we would need to make @properties *per laser*.
    # TODO: Put Getters here!

    def sanity_check(self):
        """Check if the current (live) configuration passes all pre-checks.

        It's worth calling this right before conducting an imaging run.
        """
        # Run through all checks first; raise an assertion error at the end.
        # Do the Base Class Sanity Checks first and cascade them.
        error_msgs = []
        try:
            super().sanity_check()
        except AssertionError as e:
            error_msgs.append(str(e))

        # TODO: a bunch of DISPIM-specific sanity checks.
        # TODO: put this in the base class.
        assert self.local_storage_dir.exists(), \
            f"Error: local storage directory '{self.local_storage_dir}' " \
            "does not exist."
        # Check if external storage path exists only if it was specified.
        if self.ext_storage_dir is not None:
            assert self.ext_storage_dir.exists(), \
                "Error: external storage directory " \
                f"'{self.ext_storage_dir}' does not exist."

        # Create a big error message at the end.
        if len(error_msgs):
            all_msgs = "\n".join(error_msgs)
            raise AssertionError(all_msgs)

