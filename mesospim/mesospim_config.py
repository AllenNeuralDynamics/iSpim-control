#!/usr/bin/env python3
"""Dynamic Mesospim TOML config that can be loaded, reloaded, and saved."""

import logging
import copy
from functools import cached_property
from pathlib import Path
from .config_base import TomlConfig

# Constants:
MAX_OPEN_LOOP_Z_DIST_UM = 1
# Constants below are used for sanity checks. Do not change these.
MIN_LINE_INTERVAL = 10e-6  # seconds. min time a row is turned on.
MAX_ADJ_GALVO_VOLTAGE = 10  # max differential voltage abs(a-b) that the
                            # galvo can withstand.
# TODO: ETL max voltage constants.

# A template from which we can generate a blank mesospim config toml file.
TomlTemplate = \
    {
        "imaging_params":
            {
                "local_storage_directory": ".",
                "external_storage_directory": ".",
                "subject_id": "brain_00",
                "tile_prefix": "tile",
                "tile_overlap_x_percent": 15,
                "tile_overlap_y_percent": 15,
                "z_step_size_um": 1,
                "volume_x_um": 2048,
                "volume_y_um": 2048,
                "volume_z_um": 1,
            }
        # TODO: populate remaining fields.
    }


class MesospimConfig(TomlConfig):
    """A Mesospim Configuration."""

    def __init__(self, toml_filepath: str):
        """Read config file. Warn if not found, but create sensible defaults."""
        super().__init__(toml_filepath, TomlTemplate)
        self._load_attributes()
        # Run sanity checks.
        self._sanity_check_config()

    def reload(self):
        super().reload()
        # We must call this every time we reload a cfg file.
        self._load_attributes()
        # Clear cached properties since they need to be recomputed.
        cls = self.__class__
        cached_properties = [a for a in dir(self) if
                             isinstance(getattr(cls, a, cls), cached_property)]
        for prop in cached_properties:
            # Check if hasattr first; otherwise delattr will fail on any
            # cached_property that was never cached.
            if hasattr(self, prop):
                delattr(self, prop)

    def _load_attributes(self):
        """Associate config values to attributes.

        Assigning attributes this way makes them easier to access externally
        rather than needing to know the config structure.

        Note that most of these attributes are immutable. Each time we
        load the config, we must call this function again to reassign them.

        Note that these attributes are READ-ONLY. Making them writeable means
        we would need to do an @<attribute_name>.setter for each one of them.
        """
        # NOTE: float conversion in various places is to avoid tomlkit bug
        # involving <tomlkit.items.Integer> + <float>

        # Expose these to top level since they are kwds for object __init__.
        self.etl_obj_params = self.cfg['etl_driver_kwds']
        self.tiger_obj_params = self.cfg['tiger_controller_driver_kwds']
        self.daq_obj_params = self.cfg['daq_driver_kwds']
        self.filter_wheel_params = self.cfg['filter_wheel_kwds']

        dcam_specs = self.cfg['dcam_specs']
        # Time between moving the slit by one row. (AKA: "line interval")
        self.row_interval = dcam_specs['row_interval']
        # Lightsheet scan direction: forward or backward.
        self.scan_direction = dcam_specs['scan_direction']

        # EXPERIMENT PARAMETERS loaded from file.
        imaging_specs = self.cfg['imaging_params']
        self.volume_x_um = float(imaging_specs['volume_x_um'])
        self.volume_y_um = float(imaging_specs['volume_y_um'])
        self.volume_z_um = float(imaging_specs['volume_z_um'])
        self.tile_prefix = imaging_specs['tile_prefix']  # specimen names

        design_specs = self.cfg['design_specs']
        self.tile_size_um = float(design_specs['tile_size_um'])
        self.slit_width = design_specs['slit_width']
        self.sensor_row_count = design_specs['sensor_row_count']

        self.tile_overlap_x_percent = imaging_specs['tile_overlap_x_percent']
        self.tile_overlap_y_percent = imaging_specs['tile_overlap_y_percent']
        self.subject_id = imaging_specs['subject_id']

        stage_params = self.cfg['sample_stage_specs']
        self.stage_move_time_um = stage_params['move_time_um']
        self.stage_backlash_reset_dist_um = \
            stage_params['backlash_reset_distance_um']

        waveform_specs = self.cfg['waveform_specs']
        tuning_galvo_specs = waveform_specs['adjustment_galvo']

        self.etl_duty = 1.0  # pure sawtooth

        # Waveform Specs. Pull most from the config. Derive the remaining ones.
        # additional delay for the etl
        self.start_of_frame_delay = waveform_specs['start_of_frame_delay']
        self.adj_galvo_setpoint = tuning_galvo_specs['voltage_setpoint']

        # Daq Hardware Specs
        self.daq_update_freq = self.daq_obj_params['update_frequency']

        # Laser Channel Specs
        self.laser_specs = self.cfg['channel_specs']

        # Estimates
        self.tiles_per_second = \
            float(self.cfg['estimates']['tiles_per_second'])
        self.gigabytes_per_image = \
            float(self.cfg['estimates']['gigabytes_per_image'])

    # Getters. These values depend on the laser wavelength.

    def get_etl_settling_time(self, wavelength: int):
        """Return the laser-wavelength-specific settling time for the etl."""
        # wavelength key is read as a string in the toml file.
        return self.laser_specs[str(wavelength)]['etl']['settling_time']

    def get_etl_delay(self, wavelength: int):
        # wavelength key is read as a string in the toml file.
        return self.laser_specs[str(wavelength)]['etl']['delay']

    def get_waveform_cycle_time(self, wavelength: int):
        # cycle time of waveforms not including delays
        return self.total_exposure_time + \
               self.get_etl_settling_time(wavelength)

    def get_waveform_delay(self, wavelength: int):
        # delay imposed on all waveforms (except etl) before the signals play.
        return self.get_etl_delay(wavelength) + self.start_of_frame_delay

    def get_daq_cycle_time(self, wavelength: int):
        # Total time the daq spends playing one period of the waveform.
        return self.get_waveform_cycle_time(wavelength) + \
               self.get_waveform_delay(wavelength)

    def get_laser_duty_cycle(self, wavelength: int):
        return self.total_exposure_time / \
               self.get_waveform_cycle_time(wavelength)

    def get_daq_num_samples(self, wavelength: int):
        """Total samples per stored daq channel for waveform generation."""
        return round(self.daq_update_freq*self.get_daq_cycle_time(wavelength))

    # Any derived parameter not explicitly in the config is an @property
    # so we can reload the config without needing to recompute properties.

    @property
    def imaging_wavelengths(self):
        """Returns the list of wavelengths used just for imaging in order.

        Note: this may be a subset of all configured wavelengths.
        Note: repeats are allowed since this list is interpretted as an
            execution order, but it is rare.
        """
        return self.cfg['imaging_params']['laser_wavelengths']

    @cached_property
    def xy_voxel_size_um(self):
        return float(self.tile_size_um) / self.sensor_row_count

    @cached_property
    def z_voxel_size_um(self):
        return self.cfg['imaging_params'].get('z_step_size_um',
                                              self.xy_voxel_size_um)

    @cached_property
    def local_storage_dir(self):
        return Path(self.cfg['imaging_params']['local_storage_directory'])

    @cached_property
    def ext_storage_dir(self):
        # External storage is optional and allowed to be unspecified.
        # If unspecified, create directory in the local storage directory.
        return Path(self.cfg['imaging_params'].get(
                                                'external_storage_directory',
                                                self.local_storage_dir))

    @cached_property
    def row_exposure_time(self):
        # Total time any row gets exposed to laser.
        return self.row_interval*self.slit_width

    @cached_property
    def total_exposure_time(self):
        # Total time spent where any part of the sensor is being exposed.
        return self.sensor_row_count*self.row_interval + self.row_exposure_time

    @property
    def laser_wavelengths(self):
        """Returns set of all configured laser wavelengths.

        Note: this is NOT the subset of wavelengths used for imaging."""
        return set([int(nm) for nm in self.cfg['channel_specs'].keys()])

    @property
    def daq_used_channels(self):
        """Return the total channels used on the daq."""
        # ao channels for lasers must be tallied up from channel_specs.
        # add 1 for the digital "output trigger" signal.
        ao_laser_count = 0
        # Since it's possible that lasers aren't strictly driven by an
        # ao channel, we must tally them up.
        for wavelen, specs in self.laser_specs.items():
            if 'ao_channel' in specs:
                ao_laser_count += 1
        return len(self.cfg['daq_ao_names_to_channels']) + \
            ao_laser_count + 1

    @property
    def daq_ao_names_to_channels(self):
        """Return a dict of {<analog output signal name> : <daq ao channel>}"""
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

    def _sanity_check_config(self):

        # Run through all checks first, but raise the error at the end.
        # TODO: sanity check ETL voltage limits.
        error_msgs = []

        # Check if there is at least one laser wavelength to image with.
        if len(self.imaging_wavelengths) < 1:
            msg = "At least one laser must be specified to image with."
            self.log.error(msg)
            error_msgs.append(msg)
        # Warn if there are repeat values in imaging wavelengths.
        if len(set(self.imaging_wavelengths)) > len(self.imaging_wavelengths):
            self.log.warning("Repeat values are present in the sequence of "
                             "lasers to image with.")
        # Check if z voxel size was user-specified as different from isotropic.
        z_voxel_size_um = self.cfg['imaging_params'].get('z_step_size_um', None)
        if z_voxel_size_um is None:
            self.log.warning("Z Step Size not listed. Defaulting to an "
                             f"isotropic step size of {self.xy_voxel_size_um}")
        # Warn on no external storage dir.
        if self.ext_storage_dir == self.local_storage_dir:
            self.log.warning("Output storage directory unspecified. Data will "
                             f"be saved to {self.ext_storage_dir.resolve()}.")
        # TODO: Throw error on out-of-bounds stage x, y, or z movement.
        # TODO: finish refactor.
        assert 0 < self.tile_overlap_x_percent < 100, \
            f"Error: Specified x overlap ({self.tile_overlap_x_percent}) " \
            "is out of bounds."
        assert 0 < self.tile_overlap_y_percent < 100, \
            f"Error: Specified x overlap ({self.tile_overlap_x_percent}) " \
            "is out of bounds."
        assert self.row_interval >= MIN_LINE_INTERVAL, \
            f"Error: row interval ({self.row_interval*1e6:.3f} [us]) too " \
            f"small. Minimum is {MIN_LINE_INTERVAL*1e6:.3f} [us]."
        assert abs(self.adj_galvo_setpoint * 2) < MAX_ADJ_GALVO_VOLTAGE, \
            f"Error: adjustment galvo setpoint voltage is out of range."
        assert self.scan_direction in {'FORWARD', 'BACKWARD'}, \
            "Error. Typo in config.toml scan direction."
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
            all_msgs = "\n".join(msg for msg in error_msgs)
            raise RuntimeError(all_msgs)

    def print_summary_stats(self):
        # Print some stats:
        print("--Config Stats--")
        print("  Total \"sensor-active\" time (computed): "
              f"{self.total_exposure_time*1e3:.3f} [ms]")
        #print(f"  DAQ cycle time (computed): {self.daq_cycle_time*1e3:.3f} [ms]")
        #print()
        print("Volume Capture Stat:")
        print(f"  percent overlap x: {self.tile_overlap_x_percent:.1f}%")
        print(f"  percent overlap y: {self.tile_overlap_x_percent:.1f}%")
        print("  Desired dimensions: "
              f"{self.volume_x_um:.1f}[um] x "
              f"{self.volume_y_um:.1f}[um] x "
              f"{self.volume_z_um:.1f}[um].")
        print(f"  voxel (x,y,z) size: {self.xy_voxel_size_um:.3f}[um] x "
              f"{self.xy_voxel_size_um:.3f}[um] x "
              f"{self.z_voxel_size_um:.3f}[um].")
        print()
        print("GUI settings for debugging:")
        print(f"  Line Interval: {self.row_interval*1e6:.3f} [us]")
        print(f"  Exposure Time (computed): {self.row_exposure_time*1e3:.3f} [ms]")
        print()
