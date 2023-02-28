"""Mesospim Config object to flatten item access in the TOML-based config."""

from spim_core.config_base import SpimConfig
from ispim.config_template import TomlTemplate
import copy


class IspimConfig(SpimConfig):
    """A Dispim Configuration."""

    def __init__(self, toml_filepath: str, create: bool = False):
        """Read config file. Warn if not found, but create sensible defaults."""
        super().__init__(toml_filepath, TomlTemplate, create=create)

        # Note: these are mutable, so reloading the toml doesn't affect them.
        self.imaging_specs = self.cfg['imaging_specs']
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
        # self.filter_wheel_kwds = self.cfg['filter_wheel_kwds']

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
        return round(self.daq_update_freq * self.delay_time)

    def get_rest_samples(self):
        """Return the rest samples between the left and right views."""
        return round(self.daq_update_freq * self.rest_time)

    def get_exposure_samples(self):
        """Return the exposure samples between the left and right views."""
        return round(self.daq_update_freq * self.exposure_time)

    def get_period_samples(self):
        """Return the exposure samples between the left and right views."""
        return round(self.daq_update_freq * self.get_period_time())

    def get_daq_cycle_samples(self):
        """Return the total waveform cycle time for a frame."""
        return round(self.daq_update_freq * self.get_daq_cycle_time())

    # TODO: consider putting this in the base class since literally every
    #   machine has a sample.
    @property
    def sample_pose_kwds(self):
        return self.cfg['sample_pose_kwds']

    @property
    def scan_direction(self):
        """Lightsheet scan direction: forward or backward."""
        return self.camera_specs['scan_direction']

    @scan_direction.setter
    def scan_direction_left(self, dir:str):
        # Lightsheet scan direction: forward or backward.
        self.camera_specs['scan_direction'] = dir

    @property
    def line_time(self):
        """Defines line rate of camera in us"""
        return self.waveform_specs['line_time_us']

    @line_time.setter
    def line_time(self, us: float):
        """Sets line rate of camera in us"""
        self.waveform_specs['line_time_us'] = us

    @property
    def slit_width_pix(self):
        """Returns the slit width in pixels.
        :unit px"""
        return self.design_specs['slit_width_pixels']

    @slit_width_pix.setter
    def slit_width_pix(self, width: int):
        """Sets the slit width in pixels."""
        self.design_specs['slit_width_pixels'] = width

    @property
    def sensor_row_count(self):
        return self.tile_specs['row_count_pixels']

    @property
    def sensor_column_count(self):
        """Pixels in column?"""
        return self.tile_specs['column_count_pixels']

    @property
    def daq_update_freq(self):
        """Frequency DAQ updates
        :unit hz"""
        return self.daq_obj_kwds['update_frequency_hz']

    @daq_update_freq.setter
    def daq_update_freq(self, hz: int):
        self.daq_obj_kwds['update_frequency_hz'] = hz

    # TODO: consider putting this in the parent class.
    # TODO: handle case if we want this to default to something else.

    @property
    def imaging_wavelengths(self):
        """laser wavelengths used in imaging
        :unit nm"""
        return self.cfg['imaging_specs']['laser_wavelengths']

    @imaging_wavelengths.setter
    def imaging_wavelengths(self, wl: int):
        self.cfg['imaging_specs']['laser_wavelengths'] = wl

    @property
    def z_step_size_um(self):
        """z step size in um
        :unit um"""
        return self.cfg['imaging_specs']['z_step_size_um']

    @z_step_size_um.setter
    def z_step_size_um(self, um: float):
        self.cfg['imaging_specs']['z_step_size_um'] = um

    @property
    def scan_speed_mm_s(self):
        """Return the volumetric scan speed of the stage."""
        jitter_time_s = 0.01  # 10 ms jitter time for stage pulses
        step_size_mm = self.imaging_specs['z_step_size_um'] / 1000.0
        scan_speed_mm_s = step_size_mm / (self.get_daq_cycle_time() + jitter_time_s)
        return scan_speed_mm_s

    # TODO: consider putting this in the parent class.
    @property
    def stage_backlash_reset_dist_um(self):
        """stage backlash reset distance um
        :unit um"""
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
        """Return the delay time between the left and right views.
        :unit s"""
        return self.waveform_specs['delay_time']

    @delay_time.setter
    def delay_time(self, delay_time: float):
        self.waveform_specs['delay_time'] = delay_time

    @property
    def rest_time(self):
        """Return the delay time between the left and right views.
        :unit s"""
        return self.waveform_specs['rest_time']

    @rest_time.setter
    def rest_time(self, rest_time: float):
        self.waveform_specs['rest_time'] = rest_time

    @property
    def exposure_time(self):
        """Return the total exposure time for a frame.
        :unit s"""
        return self.waveform_specs['exposure_time']

    @exposure_time.setter
    def exposure_time(self, exposure_time: float):
        self.waveform_specs['exposure_time'] = exposure_time

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
