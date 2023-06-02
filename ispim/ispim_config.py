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
        self.experiment_specs = self.cfg['experiment_specs']
        self.stage_specs = self.cfg['stage_specs']
        self.tiger_specs = self.cfg['tiger_specs']
        self.laser_specs = self.cfg['channel_specs']
        self.design_specs = self.cfg['design_specs']
        self.camera_specs = self.cfg['camera_specs']
        self.waveform_specs = self.cfg['waveform_specs']
        self.tiger_obj_kwds = self.cfg['tiger_controller_driver_kwds']
        self.daq_obj_kwds = self.cfg['daq_driver_kwds']
        self.filter_wheel_kwds = self.cfg['filter_wheel_kwds']

    # Getters. These must be functions since values depend on the laser
    # wavelength. Otherwise, we would need to make @properties *per laser*.

    def get_period_time(self):
        """Return the total waveform period time for a frame."""

        pre_buffer = self.pre_buffer_time_s if self.pre_buffer_time_s > self.laser_pre_buffer_time_s else self.laser_pre_buffer_time_s
        post_buffer = self.post_buffer_time_s if self.post_buffer_time_s > self.laser_post_buffer_time_s else self.laser_post_buffer_time_s
        return pre_buffer+ self.exposure_time + post_buffer + self.rest_time

    def get_pre_buffer_samples(self):
        """Return the buffer samples before actuating the galvos."""
        return round(self.daq_update_freq * self.pre_buffer_time_s)

    def get_post_buffer_samples(self):
        """Return the buffer samples after actuating the galvos."""
        return round(self.daq_update_freq * self.post_buffer_time_s)

    def get_rest_samples(self):
        """Return the rest samples between interleaved channels."""
        return round(self.daq_update_freq * self.rest_time)

    def get_exposure_samples(self):
        """Return the exposure samples between the left and right views."""
        return round(self.daq_update_freq * self.exposure_time)

    def get_period_samples(self):
        """Return the exposure samples between the left and right views."""
        return round(self.daq_update_freq * self.get_period_time())

    def get_laser_post_buffer_samples(self):
        """Return the buffer samples of laser delay."""
        return round(self.daq_update_freq * self.laser_post_buffer_time_s)

    def get_laser_pre_buffer_samples(self):
        """Return the buffer samples of laser delay."""
        return round(self.daq_update_freq * self.laser_pre_buffer_time_s)

    @property
    def acquisition_style(self):
        """Returns whether acquisition will play interleaved waveforms at each tile
        or take sequential rounds of tiling"""
        # TODO: Error if not interleaved or sequential?
        return self.imaging_specs['acquisition_style']

    # TODO: consider putting this in the base class since literally every
    #   machine has a sample.
    @property
    def sample_pose_kwds(self):
        return self.cfg['sample_pose_kwds']

    @property
    def experimenters_name(self):
        return self.experiment_specs['experimenters_name']

    @experimenters_name.setter
    def experimenters_name(self, name: str):
        self.experiment_specs['experimenters_name'] = name

    @property
    def immersion_medium(self):
        return self.experiment_specs['immersion_medium']

    @immersion_medium.setter
    def immersion_medium(self, medium : str):
        self.experiment_specs['immersion_medium'] = medium

    @property
    def immersion_medium_refractive_index(self):
        return self.experiment_specs['immersion_medium_refractive_index']

    @immersion_medium_refractive_index.setter
    def immersion_medium_ri(self, ri: float):
        self.experiment_specs['immersion_medium_refractive_index'] = ri

    @property
    def scan_direction(self):
        """Lightsheet scan direction: forward or backward."""
        return self.camera_specs['scan_direction']

    @scan_direction.setter
    def scan_direction(self, direction: str):
        """Sets line rate of camera in us"""
        self.camera_specs['scan_direction'] = direction

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
    def jitter_time_s(self):
        """jitter time for stage
        :unit s"""
        return self.stage_specs['jitter_time_s']

    @jitter_time_s.setter
    def jitter_speed_s(self, seconds: float):
        self.stage_specs['jitter_time_s'] = seconds

    @property
    def scan_speed_mm_s(self):
        """Return the volumetric scan speed of the stage."""

        step_size_mm = self.imaging_specs['z_step_size_um'] / 1000.0
        if self.acquisition_style == 'interleaved':
            scan_speed_mm_s = (
                        step_size_mm / ((self.get_period_time() * len(self.imaging_wavelengths)) + self.jitter_time_s))
        elif self.acquisition_style == 'sequential':
            scan_speed_mm_s = (step_size_mm / ((self.get_period_time()) + self.jitter_time_s))
        # TODO: Error if niether one?
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
    def post_buffer_time_s(self):
        """Return the buffer time after actuating the galvos.
        :unit s"""
        return self.waveform_specs['post_buffer_time_s']

    @post_buffer_time_s.setter
    def post_buffer_time_s(self, post_buffer_time_s: float):
        self.waveform_specs['post_buffer_time_s'] = post_buffer_time_s

    @property
    def pre_buffer_time_s(self):
        """Return the buffer time before actuating the galvos.
        :unit s"""
        return self.waveform_specs['pre_buffer_time_s']

    @pre_buffer_time_s.setter
    def pre_buffer_time_s(self, pre_buffer_time_s: float):
        self.waveform_specs['pre_buffer_time_s'] = pre_buffer_time_s

    @property
    def laser_pre_buffer_time_s(self):
        """Return the buffer time before actuating the galvos.
        :unit s"""
        return self.waveform_specs['laser_pre_buffer_time_s']

    @laser_pre_buffer_time_s.setter
    def laser_pre_buffer_time_s(self, pre_buffer_time_s: float):
        self.waveform_specs['laser_pre_buffer_time_s'] = pre_buffer_time_s

    @property
    def laser_post_buffer_time_s(self):
        """Return the buffer time before actuating the galvos.
        :unit s"""
        return self.waveform_specs['laser_post_buffer_time_s']

    @laser_post_buffer_time_s.setter
    def laser_post_buffer_time_s(self, post_buffer_time_s: float):
        self.waveform_specs['laser_post_buffer_time_s'] = post_buffer_time_s

    @property
    def rest_time(self):
        """Return the rest time between interleaved channels.
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
        return set([int(nm) for nm in self.cfg['channel_specs'].keys() if nm.isdigit()])

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

