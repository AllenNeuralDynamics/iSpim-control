"""Mesospim Config object to flatten item access in the TOML-based config."""

from mesospim.config_base import SpimConfig
import copy

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
                "z_step_size_um": 1,
                "volume_x_um": 2304,
                "volume_y_um": 2304,
                "volume_z_um": 1,
            }
        # TODO: populate remaining fields.
    }


class DispimConfig(SpimConfig):
    """A Dispim Configuration."""

    def __init__(self, toml_filepath: str):
        """Read config file. Warn if not found, but create sensible defaults."""
        super().__init__(toml_filepath, TomlTemplate)

        # Note: these are mutable, so reloading the toml doesn't affect them.
        self.stage_specs = self.cfg['stage_specs']
        self.tiger_specs = self.cfg['tiger_specs']
        self.laser_specs = self.cfg['channel_specs']
        self.design_specs = self.cfg['design_specs']
        self.camera_specs = self.cfg['camera_specs']
        self.waveform_specs = self.cfg['waveform_specs']
        self.tiger_obj_kwds = self.cfg['tiger_controller_driver_kwds']
        self.daq_obj_kwds = self.cfg['daq_driver_kwds']
        # TODO: dispim has 2 filterwheels. We must set the location of both
        #   programmatically.
        #self.filter_wheel_kwds = self.cfg['filter_wheel_kwds']

    # Getters. These must be functions since values depend on the laser
    # wavelength. Otherwise, we would need to make @properties *per laser*.

    def get_period_time(self):
        """Return the total waveform cycle time for a frame."""
        return self.exposure_time + self.rest_time

    def get_daq_cycle_time(self):
        """Return the total waveform cycle time for a frame."""
        return self.get_period_time() + self.delay_time

    def get_delay_samples(self):
        """Return the delay samples between the left and right views."""
        return round(self.daq_update_freq*self.delay_time)

    def get_rest_samples(self):
        """Return the rest samples between the left and right views."""
        return round(self.daq_update_freq*self.rest_time)

    def get_exposure_samples(self):
        """Return the exposure samples between the left and right views."""
        return round(self.daq_update_freq*self.exposure_time)

    def get_period_samples(self):
        """Return the exposure samples between the left and right views."""
        return round(self.daq_update_freq*self.get_period_time())

    def get_daq_cycle_samples(self):
        """Return the total waveform cycle time for a frame."""
        return round(self.daq_update_freq*self.get_daq_cycle_time())

    # TODO: consider putting this in the base class since literally every
    #   machine has a sample.
    @property
    def sample_pose_kwds(self):
        return self.cfg['sample_pose_kwds']

    # @property
    # def row_interval(self):
    #     # Time between moving the slit by one row. (AKA: "line interval")
    #     return self.camera_specs['row_interval']

    # @property
    # def scan_direction_left(self):
    #     # Lightsheet scan direction: forward or backward.
    #     return self.camera_specs['camera_left']['scan_direction']

    # @property
    # def scan_direction_left(self):
    #     # Lightsheet scan direction: forward or backward.
    #     return self.camera_specs['camera_right']['scan_direction']

    # @scan_direction_left.setter
    # def scan_direction(self, direction: DCAMPROP.READOUT_DIRECTION):
    #     self.dcam_specs['scan_direction'] = direction.name

    # @property
    # def slit_width(self):
    #     """Returns the slit width in pixels."""
    #     return self.design_specs['slit_width']

    # @slit_width.setter
    # def slit_width(self, width: int):
    #     """Sets the slit width in pixels."""
    #     self.design_specs['slit_width'] = width

    @property
    def sensor_row_count(self):
        return self.tile_specs['row_count_pixels']

    @property
    def sensor_column_count(self):
        return self.tile_specs['column_count_pixels']

    # @property
    # def start_of_frame_delay(self):
    #     return self.cfg['waveform_specs']['start_of_frame_delay']

    # @start_of_frame_delay.setter
    # def start_of_frame_delay(self, seconds: float):
    #     self.cfg['waveform_specs']['start_of_frame_delay'] = seconds

    # @property
    # def trim_galvo_setpoint(self):
    #     return self.cfg['waveform_specs']['trim_galvo']['voltage_setpoint']

    # @trim_galvo_setpoint.setter
    # def trim_galvo_setpoint(self, volts: float):
    #     self.cfg['waveform_specs']['trim_galvo']['voltage_setpoint'] = volts

    @property
    def daq_update_freq(self):
        return self.daq_obj_kwds['update_frequency_hz']

    @daq_update_freq.setter
    def daq_update_freq(self, hz: int):
        self.daq_obj_kwds['update_frequency_hz'] = hz

    # TODO: consider putting this in the parent class.
    # TODO: handle case if we want this to default to something else.
    @property
    def z_step_size_um(self):
        return self.cfg['imaging_specs']['z_step_size_um']

    @z_step_size_um.setter
    def z_step_size_um(self, um: float):
        self.cfg['imaging_specs']['z_step_size_um'] = um

    # TODO: consider putting this in the parent class.
    @property
    def stage_backlash_reset_dist_um(self):
        return self.stage_specs['backlash_reset_distance_um']

    @stage_backlash_reset_dist_um.setter
    def stage_backlash_reset_dist_um(self, micrometers: int):
        self.stage_specs['backlash_reset_distance_um'] = micrometers

    # Any derived parameter not explicitly in the config is an @property,
    # so we can reload the config without needing to recompute properties.
    # These do NOT get setters.
    # FIXME: make separate getters for XY.
    #   Handle Overlap correctly.

    @property
    def delay_time(self):
        """Return the delay time between the left and right views."""
        return self.waveform_specs['delay_time']

    @delay_time.setter
    def delay_time(self, delay_time: float):
        self.waveform_specs['delay_time'] = delay_time

    @property
    def rest_time(self):
        """Return the delay time between the left and right views."""
        return self.waveform_specs['rest_time']

    @rest_time.setter
    def rest_time(self, rest_time: float):
        self.waveform_specs['rest_time'] = rest_time

    @property
    def exposure_time(self):
        """Return the total exposure time for a frame."""
        return self.waveform_specs['exposure_time']

    @exposure_time.setter
    def exposure_time(self, exposure_time: float):
        self.waveform_specs['exposure_time'] = exposure_time

    @property
    def camera_left_offset(self):
        return self.camera_specs['camera_left']['offset']

    @camera_left_offset.setter
    def camera_left_offset(self, left_offset: float ):
        self.camera_specs['camera_left']['offset'] = left_offset

    @property
    def camera_right_offset(self):
        return self.camera_specs['camera_right']['offset']

    @camera_right_offset.setter
    def camera_right_offset(self, right_offset: float):
        self.camera_specs['camera_right']['offset'] = right_offset

    @property
    def laser_wavelengths(self):
        """Returns set of all configured laser wavelengths.
        Note: this is NOT the subset of wavelengths used for imaging."""
        return set([int(nm) for nm in self.cfg['channel_specs'].keys()])

    @property
    def daq_used_channels(self):
        """Return the total channels used on the daq."""
        # ao channels for lasers must be tallied up from channel_specs.
        ao_laser_count = 0
        # Since it's possible that lasers aren't strictly driven by an
        # ao channel, we must tally them up.
        for wavelength, specs in self.laser_specs.items():
            if 'ao_channel' in specs:
                ao_laser_count += 1
        return len(self.cfg['daq_ao_names_to_channels']) + ao_laser_count

    @property
    def daq_ao_names_to_channels(self):
        """Return a dict of {<analog output signal name> : <daq ao channel>}.

        Laser signals are stuffed in as str(wavelength_in_nm).
        """
        # Since this data is generated from mutable values, make a deepcopy
        # so we don't change TOML values if we later save the TOML.
        ao_names_to_channels = \
            copy.deepcopy(self.cfg['daq_ao_names_to_channels'])
        for wavelen, specs in self.laser_specs.items():
            ao_channel = specs.get('ao_channel', None)
            # Handle (rare!) case that a laser isn't driven by an ao channel.
            if ao_channel is not None:
                ao_names_to_channels[f"{wavelen}"] = ao_channel
            # This case is so rare that we should warn about it.
            else:
                self.log.warning(f"{wavelen} [nm] laser is not driven by an"
                                 f"analog output channel on the NI DAQ.")
        return ao_names_to_channels

    # Simple @properties that do not have setters because they shouldn't be
    # changed from anything else other than the toml file itself.
    @property
    def tiles_per_second(self):
        return float(self.cfg['estimates']['tiles_per_second'])

    # DO WE NEED THIS?
    # def sanity_check(self):
    #     """Check if the current (live) configuration passes all pre-checks.
    #     It's worth calling this right before conducting an imaging run.
    #     """
    #     # Run through all checks first; raise an assertion error at the end.
    #     # Do the Base Class Sanity Checks first and cascade them.
    #     error_msgs = []
    #     try:
    #         super().sanity_check()
    #     except AssertionError as e:
    #         error_msgs.append(str(e))

    #     # TODO: a bunch of DISPIM-specific sanity checks.
    #     # TODO: put this in the base class.
    #     assert self.local_storage_dir.exists(), \
    #         f"Error: local storage directory '{self.local_storage_dir}' " \
    #         "does not exist."
    #     # Check if external storage path exists only if it was specified.
    #     if self.ext_storage_dir is not None:
    #         assert self.ext_storage_dir.exists(), \
    #             "Error: external storage directory " \
    #             f"'{self.ext_storage_dir}' does not exist."

    #     # Create a big error message at the end.
    #     if len(error_msgs):
    #         all_msgs = "\n".join(error_msgs)
    #         raise AssertionError(all_msgs)

    def print_summary_stats(self):
        # Print some stats:
        print("--Config Stats--") # TO DO ADD DISPIM VALUES HERE
        print("Volume Capture Stat:")
        print(f"  percent overlap x: {self.tile_overlap_x_percent:.1f}%")
        print(f"  percent overlap y: {self.tile_overlap_x_percent:.1f}%")
        print("  Desired dimensions: "
              f"{self.volume_x_um:.1f}[um] x "
              f"{self.volume_y_um:.1f}[um] x "
              f"{self.volume_z_um:.1f}[um].")
        print()
        print("GUI settings for debugging:") # TO DO ADD DISPIM VALUES HERE
        print()