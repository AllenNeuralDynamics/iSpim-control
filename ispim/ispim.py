"""Abstraction of the ispim Instrument."""
from datetime import timedelta, datetime
import calendar
import logging
import numpy as np
from pathlib import Path
from time import perf_counter, sleep, time
from mock import NonCallableMock as Mock
from threading import Thread, Event
from ispim.ispim_config import IspimConfig
from ispim.devices.ni import WaveformHardware
from ispim.compute_waveforms import generate_waveforms
from tigerasi.tiger_controller import TigerController, STEPS_PER_UM
from tigerasi.device_codes import PiezoControlMode, TTLIn0Mode
from tigerasi.sim_tiger_controller import SimTigerController as SimTiger
from spim_core.spim_base import Spim
from spim_core.devices.tiger_components import SamplePose,FilterWheel
from spim_core.processes.data_transfer import DataTransfer
import os
from acquire import DeviceState
import tifffile
import shutil
#from vortran_laser import stradus
from ispim.operations import normalized_dct_shannon_entropy
import threading
import sys
class Ispim(Spim):

    def __init__(self, config_filepath: str,
                 simulated: bool = False):

        self.log = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        # Log setup is handled in the parent class if we pass in a logger.
        super().__init__(config_filepath, simulated=simulated)

        self.cfg = IspimConfig(config_filepath)
        # Instantiate hardware devices
        self.ni = WaveformHardware(**self.cfg.daq_obj_kwds) if not self.simulated else \
            Mock(WaveformHardware)
        self.tigerbox = TigerController(**self.cfg.tiger_obj_kwds) if not \
            self.simulated else SimTiger(**self.cfg.tiger_obj_kwds,
                                         build_config={'Motor Axes': ['X', 'Y', 'Z', 'V', 'W', 'A', 'B', 'C', 'D']})
        self.sample_pose = SamplePose(self.tigerbox, **self.cfg.sample_pose_kwds)
        self.filter_wheel = FilterWheel(self.tigerbox,
                                        **self.cfg.filter_wheel_kwds)
        self.camera = None
        self.lasers = {}  # populated in _setup_lasers.
        self.channel_gene = {}  # dictionary containing labeled gene for each channel

        # Extra Internal State attributes for the current image capture
        # sequence. These really only need to persist for logging purposes.
        self.start_time = None
        self.total_tiles = None  # tiles to be captured.
        self.x_y_tiles = None    # tiles in x and y to be captured
        self.tiles_acquired = 0
        self.tile_time_s = 0
        self.est_run_time = None
        self.stage_x_pos = None
        self.stage_y_pos = None
        self.scout_mode = False

        # Setup hardware according to the config.
        self._setup_camera()
        self._setup_lasers()
        self._setup_motion_stage()
        # TODO, note NIDAQ is channel specific and gets instantiated within imaging loop

        # Internal state attributes.
        self.active_lasers = None  # Bookkeep which laser is configured.
        self.live_status = None  # Bookkeep if we are running in live mode.

        self.livestream_worker = None  # captures images during livestream
        self.livestream_enabled = Event()
        self.setting_up_livestream = False

        self.latest_frame = None
        self.latest_frame_layer = 0

        # start position of scan
        self.start_pos = None
        self.im = None

        self.stack = []
        self.overview_stack = None
        self.overview_process = None
        self.overview_set = Event()
        self.overview_imgs = []

        self.__sim_counter_count = 0

        self.stage_query_lock = threading.Lock()

    def _setup_camera(self):
        """Configure general settings and set camera settings to those specified in config"""

        self.log.debug(f"Setting up camera")
        for camera, specs in self.cfg.camera_specs.items():
            __import__(specs['driver'])
            camera_class = getattr(sys.modules[specs['driver']], specs['module'])
            kwds = dict(specs['kwds'])
            for k, v in kwds.items():
                if str(v).split('.')[0] in dir(sys.modules[specs['driver']]):
                    arg_class = getattr(sys.modules[specs['driver']], v.split('.')[0])
                    kwds[k] = getattr(arg_class, '.'.join(v.split('.')[1:]))
                else:
                    kwds[k] = eval(v) if '.' in str(v) else v
            self.camera = camera_class(**kwds) if not self.simulated else Mock()
            self.camera.configure(**specs['configure'])

    def _setup_lasers(self):
        """Setup lasers that will be used for imaging. Warm them up, etc."""

        self.log.debug(f"Setting up lasers")
        for wl, specs in self.cfg.laser_specs.items():
            if 'port' in specs['kwds'].keys() and specs['kwds']['port'] == 'COMxx':
                self.log.warning(f'Skipping setup for laser {wl} due to no COM port specified')
                continue
            __import__(specs['driver'])
            laser_class = getattr(sys.modules[specs['driver']], specs['module'])
            kwds = dict(specs['kwds'])
            for k, v in kwds.items():
                if str(v).split('.')[0] in dir(sys.modules[specs['driver']]):
                    arg_class = getattr(sys.modules[specs['driver']], v.split('.')[0])
                    kwds[k] = getattr(arg_class, '.'.join(v.split('.')[1:]))
                else:
                    kwds[k] = eval(v) if '.' in str(v) else v

            self.lasers[wl] = laser_class(**kwds) if not self.simulated else Mock()
            self.lasers[wl].disable_cdrh()  # disable five second cdrh delay
            self.log.debug(f"Successfully setup {wl} laser")

    def _setup_motion_stage(self):
        """Configure the sample stage for the ispim according to the config."""
        self.log.info("Setting backlash in Z to 0")
        self.sample_pose.set_axis_backlash(Z=0.0)
        self.log.info("Setting speeds to 1.0 mm/sec")
        self.tigerbox.set_speed(X=1.0, Y=1.0, Z=1.0)
        # Note: Tiger X is Tiling Z, Tiger Y is Tiling X, Tiger Z is Tiling Y.
        #   This axis remapping is handled upon SamplePose __init__.
        # loop over axes and verify in external mode
        # TODO, think about where to store this mapping in config
        # TODO, merge ispim commands in tigerasi
        # TODO, how to call this? via tigerbox?
        externally_controlled_axes = \
            {a.lower(): PiezoControlMode.EXTERNAL_CLOSED_LOOP for a in
             self.cfg.tiger_specs['axes'].values()}
        self.tigerbox.set_axis_control_mode(**externally_controlled_axes)

        # TODO, this needs to be buried somewhere else
        # TODO, how to store card # mappings, in config?
        (self.tigerbox
         .set_ttl_pin_modes(in0_mode=TTLIn0Mode.MOVE_TO_NEXT_ABS_POSITION,
                                        card_address=31))

    def _setup_waveform_hardware(self, active_wavelength: list, live: bool = False, scout_mode: bool = False):

        if self.simulated:
            self.ni.counter_task = Mock()
            self.ni.counter_task.read = self.__sim_counter_read
        if not self.livestream_enabled.is_set():       # Only configures daq on the initiation of livestream
            self.log.info("Configuring NIDAQ")
            self.ni.configure(self.cfg.get_period_time(), self.cfg.daq_ao_names_to_channels, len(active_wavelength), live)
        self.log.info("Generating waveforms.")
        _, voltages_t = generate_waveforms(self.cfg, active_wavelength)
        self.log.info("Writing waveforms to hardware.")
        self.ni.assign_waveforms(voltages_t, scout_mode)

    def __sim_counter_read(self):
        count = self.__sim_counter_count
        self.__sim_counter_count += 1
        return count

    # TODO: this should be a base class thing.s
    def check_ext_disk_space(self, xtiles, ytiles, ztiles):
        """Checks ext disk space before scan to see if disk has enough space scan"""
        # One tile (tiff) is ~10368 kb
        if self.cfg.imaging_specs['filetype'] == 'Tiff':
            est_stack_filesize = self.cfg.bytes_per_image * ztiles
            est_scan_filesize = est_stack_filesize*xtiles*ytiles
            if est_scan_filesize >= shutil.disk_usage(self.cfg.ext_storage_dir).free:
                self.log.error("Not enough space in external directory")
                raise
        elif self.cfg.imaging_specs['filetype'] != 'Tiff':
            self.log.warning("Checking disk space not implemented. "
                             "Proceed at your own risk")

    def check_local_disk_space(self, z_tiles):
        """Checks local disk space before scan to see if disk has enough space for two stacks"""

        #One tile (tiff) is ~10368 kb
        if self.cfg.imaging_specs['filetype'] == 'Tiff':
            est_filesize = self.cfg.bytes_per_image*z_tiles
            if est_filesize*2 >= shutil.disk_usage(self.cfg.local_storage_dir).free:
                self.log.error("Not enough space on disk. Is the recycle bin empty?")
                raise
        elif self.cfg.imaging_specs['filetype'] != 'Tiff':
            self.log.warning("Checking disk space not implemented. "
                             "Proceed at your own risk")

    def acquisition_time(self, xtiles, ytiles, ztiles):

        """Calculate acquisition time based on xtiles, ytiles, and ztiles.
        ztiles should be not interleaved """

        x_y_tiles = xtiles*ytiles
        stack_time_s = ((self.cfg.get_period_time() * len(self.cfg.imaging_wavelengths)) + self.cfg.jitter_time_s) * ztiles
        est_filesize = self.cfg.bytes_per_image * ztiles
        transfer_speed_s = self.cfg.estimates['network_speed_Bps']
        file_transfer_time_s = est_filesize/transfer_speed_s

        if file_transfer_time_s > stack_time_s and not self.overview_set.is_set():
            total_time_s = file_transfer_time_s*x_y_tiles
        else:
            total_time_s = (stack_time_s*x_y_tiles) + file_transfer_time_s
            # Add one file_transfer_time to account for last tile 
        total_time_day = total_time_s / 86400
        self.log.info(f"Scan will take approximately {total_time_day}")
        completion_date = datetime.now() + timedelta(days=total_time_day)
        date_str = completion_date.strftime("%d %b, %Y at %H:%M %p")
        weekday = calendar.day_name[completion_date.weekday()]
        self.log.info(f"Time Esimate: {total_time_day:.2f} days. Imaging run "
                      f"should finish after {weekday}, {date_str}.")
        return total_time_day

    def wait_to_stop(self, axis: str, desired_position: int):
        """Wait for stage to stop moving. IN SAMPLE POSE"""
        start = time()
        while self.sample_pose.is_moving():
            pos = self.sample_pose.get_position()
            distance = abs(pos[axis.lower()] - desired_position)
            if distance < 2.0 or time()-start > 60:
                self.tigerbox.halt()
                break
            else:
                self.log.info(f"Stage is still moving! {axis} = {pos[axis.lower()]} -> {desired_position}")
                sleep(0.5)

    def log_stack_acquisition_params(self, curr_tile_index, stack_name,
                                     z_step_size_um):
        """helper function in main acquisition loop to log the current state
        before capturing a stack of images per channel."""

        for laser in self.active_lasers:
            laser = str(laser)
            tile_schema_params = \
                {
                    'tile_number': curr_tile_index,
                    'file_name': stack_name[0],
                    'coordinate_transformations': [
                        {'scale': [self.cfg.tile_size_x_um / self.cfg.sensor_column_count,
                                                        self.cfg.tile_size_y_um / self.cfg.sensor_row_count,
                                                        z_step_size_um]},
                        {'translation':[self.stage_x_pos * 0.001,
                                                                   self.stage_y_pos * 0.001,
                                                                   self.stage_z_pos * 0.001]}
                    ],
                    'channel' : {'channel_name': self.channel_gene[laser] if laser in self.channel_gene.keys() else '',
                                     'light_source_name': self.channel_gene[laser] if laser in self.channel_gene.keys() else '',
                                             'excitation_wavelength': laser,
                                             'excitation_power': self.lasers[laser].get_setpoint(),
                                             'filter_wheel_index': '0' if self.cfg.acquisition_style == 'interleaved' else self.cfg.laser_specs[laser]["filter_index"],
                                             'filter_names': [],
                                             'detector_name' : ''
                                                },
                    'channel_name': f'{laser}',
                    'x_voxel_size': self.cfg.tile_size_x_um / self.cfg.sensor_column_count,
                    'y_voxel_size': self.cfg.tile_size_y_um / self.cfg.sensor_row_count,
                    'z_voxel_size': z_step_size_um,
                    'voxel_size_units': 'micrometers',
                    'tile_x_position': self.stage_x_pos * 0.001,
                    'tile_y_position': self.stage_y_pos * 0.001,
                    'tile_z_position': self.stage_z_pos * 0.001,
                    'tile_position_units': 'millimeters',
                    'lightsheet_angle': 45,
                    'lightsheet_angle_units': 'degrees',
                    'laser_wavelength': laser,
                    'laser_wavelength_units': "nanometers",
                    'laser_gene_name': self.channel_gene[laser] if laser in self.channel_gene.keys() else None,
                    'laser_power': self.lasers[laser].get_setpoint(),
                    'laser_power_units': 'milliwatts' if self.cfg.laser_specs[laser]['intensity_mode'] == 'power' else 'percent',
                    'filter_wheel_index': '0' if self.cfg.acquisition_style == 'interleaved' else self.cfg.laser_specs[laser]["filter_index"],
                    'tags': ['schema']
                }
            self.log.info('tile data', extra=tile_schema_params)
            settings_schema_data = \
                {'tags': ['schema']}
            # Every variable in calculate waveforms
            for key in self.cfg.laser_specs[laser]['etl']:
                settings_schema_data[f'daq_etl {key}'] = self.cfg.laser_specs[laser]["etl"][key]
            for key in self.cfg.laser_specs[laser]['galvo']:
                settings_schema_data[f'daq galvo {key}'] = self.cfg.laser_specs[laser]["galvo"][key]
            self.log.info(f'laser channel {laser} acquisition settings',
                          extra=settings_schema_data)

    def run_from_config(self):

        if self.livestream_enabled.is_set():
            self.stop_livestream()
        self.collect_volumetric_image(self.cfg.volume_x_um,
                                      self.cfg.volume_y_um,
                                      self.cfg.volume_z_um,
                                      self.cfg.z_step_size_um,
                                      self.cfg.imaging_specs['laser_wavelengths'],
                                      self.cfg.scan_speed_mm_s,
                                      self.cfg.tile_overlap_x_percent,
                                      self.cfg.tile_overlap_y_percent,
                                      self.cfg.tile_prefix,
                                      self.cfg.imaging_specs['filetype'],
                                      self.local_storage_dir,
                                      self.img_storage_dir,
                                      self.deriv_storage_dir,
                                      self.cfg.acquisition_style)

    def collect_volumetric_image(self, volume_x_um: float, volume_y_um: float,
                                 volume_z_um: float,
                                 z_step_size_um: float,
                                 channels: list,
                                 scan_speed_mm_s,
                                 tile_overlap_x_percent: float,
                                 tile_overlap_y_percent: float,
                                 tile_prefix: str,
                                 filetype: str,
                                 local_storage_dir: Path = Path("."),
                                 img_storage_dir: Path = None,
                                 deriv_storage_dir: Path = None,
                                 acquisition_style: str = 'sequential'):
        """Collect a tiled volumetric image with specified size/overlap specs.
        """

        x_grid_step_um, y_grid_step_um = self.get_xy_grid_step(tile_overlap_x_percent,
                                                               tile_overlap_y_percent)

        # Calculate number of tiles in XYZ
        # Always round up so that we cover the desired imaging region.
        xtiles, ytiles, ztiles = self.get_tile_counts(tile_overlap_x_percent,
                                                      tile_overlap_y_percent,
                                                      z_step_size_um,
                                                      volume_x_um,
                                                      volume_y_um,
                                                      volume_z_um)

        self.total_tiles = xtiles * ytiles * ztiles
        self.x_y_tiles = xtiles * ytiles
        self.log.info(f"Total tiles: {self.total_tiles}.")
        actual_vol_x_um = self.cfg.tile_size_x_um + (xtiles - 1) * x_grid_step_um
        actual_vol_y_um = self.cfg.tile_size_y_um + (ytiles - 1) * y_grid_step_um
        actual_vol_z_um = z_step_size_um * ((ztiles - 1))
        self.log.info(f"Actual dimensions: {actual_vol_x_um:.1f}[um] x "
                      f"{actual_vol_y_um:.1f}[um] x {actual_vol_z_um:.1f}[um]")
        self.log.info(f"X grid step: {x_grid_step_um} [um]")
        self.log.info(f"Y grid step: {y_grid_step_um} [um]")
        self.log.info(f"Z grid step: {z_step_size_um} [um]")
        self.log.info(f'xtiles: {xtiles}, ytiles: {ytiles}, ztiles: {ztiles}')

        # Check to see if disk has enough space for two tiles
        frames = (ztiles * len(self.cfg.imaging_wavelengths)) if acquisition_style == 'interleaved' else ztiles
        self.check_local_disk_space(frames)

        #Check if external disk has enough space
        if not self.overview_set.is_set():
            self.check_ext_disk_space(xtiles, ytiles, frames)

        # Est time scan will finish
        self.start_time = datetime.now()
        self.est_run_time = self.acquisition_time(xtiles, ytiles, ztiles)

        # Move sample to preset starting position
        if self.start_pos is not None:
            # Start position is in SAMPLE POSE
            self.log.info(f'Moving to starting position at {self.start_pos["x"]}, '
                          f'{self.start_pos["y"]}, '
                          f'{self.start_pos["z"]}')
            self.sample_pose.move_absolute(x=self.start_pos['x'], wait=False)
            self.wait_to_stop('x', self.start_pos['x'])  # wait_to_stop uses SAMPLE POSE
            self.sample_pose.move_absolute(y=self.start_pos['y'], wait=False)
            self.wait_to_stop('y', self.start_pos['y'])
            self.sample_pose.move_absolute(z=self.start_pos['z'], wait=False)
            self.wait_to_stop('z', self.start_pos['z'])
            self.log.info(f'Stage moved to {self.sample_pose.get_position()}')
        else:
            self.set_scan_start(self.sample_pose.get_position())

        transfer_processes = None  # Reference to external tiff transfer process.
        # Stage positions in SAMPLE POSE
        self.stage_x_pos, self.stage_y_pos, self.stage_z_pos = (self.start_pos['x'],
                                                                self.start_pos['y'],
                                                                self.start_pos['z'])


        # Logging for JSON schema
        acquisition_params = {'session_start_time': datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
                              'local_storage_directory': str(local_storage_dir),
                              'external_storage_directory': str(img_storage_dir),
                              'specimen_id': self.cfg.imaging_specs["subject_id"],
                              'subject_id': self.cfg.imaging_specs['subject_id'],
                              'chamber_immersion': {'medium': self.cfg.immersion_medium,
                                                                 'refractive_index': self.cfg.immersion_medium_refractive_index},
                              'instrument_id': 'iSpim 1',
                              'experimenter_full_name': [self.cfg.experimenters_name],  # Needs to be in list for AIND Schema,
                              'tags': ['schema']}
        self.log.info("acquisition parameters", extra=acquisition_params)

        try:
            for j in range(ytiles):
                # move back to x=0 which maps to z=0
                self.stage_x_pos = self.start_pos['x']  # Both in SAMPLE POSE

                self.log.info("Setting speed in Y to 1.0 mm/sec")
                self.tigerbox.set_speed(Z=1.0)  # Z maps to Y

                self.log.info(f"Moving to Y = {self.stage_y_pos}.")
                self.tigerbox.move_absolute(z=round(self.stage_y_pos), wait=False)
                self.wait_to_stop('y', self.stage_y_pos)  # Use in case stage gets stuck , wait_to_stop uses SAMPLE POSE

                for i in range(xtiles):
                    # Move to specified X position
                    self.log.debug("Setting speed in X to 1.0 mm/sec")
                    self.tigerbox.set_speed(Y=1.0)  # Y maps to X
                    self.log.debug(f"Moving to X = {round(self.stage_x_pos)}.")
                    self.tigerbox.move_absolute(y=round(self.stage_x_pos), wait=False)
                    self.wait_to_stop('x', self.stage_x_pos)  # wait_to_stop uses SAMPLE POSE

                    # If sequential, loop through k for each of active_wavelenghths and feed in list as [[wl]]
                    # e.g. [[488],[561]]. Waveform generator is expecting list so give list of one wl if sequential.
                    # If Interleaved, loop through once and send one list of all wavelengths

                    loops = 1 if acquisition_style == 'interleaved' else len(channels)
                    channel = [channels] if acquisition_style == 'interleaved' else [[wl] for wl in channels]

                    for k in range(0, loops):
                        tile_start = time()
                        # Move to specified Z position
                        self.log.debug("Setting speed in Z to 1.0 mm/sec")
                        self.tigerbox.set_speed(X=1.0)  # X maps to Z
                        self.log.debug("Applying extra move to take out backlash.")
                        z_backup_pos = -STEPS_PER_UM * self.cfg.stage_backlash_reset_dist_um
                        self.tigerbox.move_absolute(x=round(z_backup_pos))
                        self.log.info(f"Moving to Z = {self.stage_z_pos}.")
                        self.tigerbox.move_absolute(x=self.stage_z_pos, wait=False)
                        self.wait_to_stop('z', self.stage_z_pos)  # wait_to_stop uses SAMPLE POSE

                        self.log.info(f"Setting scan speed in Z to {scan_speed_mm_s} mm/sec.")
                        self.tigerbox.set_speed(X=scan_speed_mm_s)
                        self.log.info(f"Actual speed {self.tigerbox.get_speed('x')}mm/sec.")

                        self.log.info(f"Setting up lasers for active channels: {channel[k]}")
                        self.setup_imaging_for_laser(channel[k])

                        # Setup capture of next Z stack.

                        # if filetype is trash for overview, it'll be zarr but doesn't matter
                        filetype_suffix = 'tiff' if filetype == 'Tiff' else 'zarr'

                        channel_string = '_'.join(map(str, self.active_lasers)) if \
                            acquisition_style == 'interleaved' else channel[k][0]
                        filenames = [f"{tile_prefix}_X_{i:0>4d}_Y_{j:0>4d}_Z_{0:0>4d}_ch_{channel_string}.{filetype_suffix}"]


                        os.makedirs(local_storage_dir, exist_ok=True)  # Make local directory if not already created
                        filepath_srcs = [local_storage_dir / f for f in filenames]

                        self.log.info(f"Collecting tile stacks at "
                                      f"({self.stage_x_pos / STEPS_PER_UM}, "
                                      f"{self.stage_y_pos / STEPS_PER_UM}) [um] "
                                      f"for channels {channel[k]} and saving to: {filepath_srcs}")
                        self.log_stack_acquisition_params(self.tiles_acquired,
                                                          filenames,
                                                          z_step_size_um)
                        # Convert to [mm] units for tigerbox.
                        slow_scan_axis_position = self.stage_x_pos / STEPS_PER_UM / 1000.0
                        self._collect_stacked_tiff(slow_scan_axis_position,
                                                   ztiles,
                                                   z_step_size_um,
                                                   filepath_srcs,
                                                   filetype,
                                                   acquisition_style)

                        # Start transferring file to its destination.
                        # Note: Image transfer is faster than image capture, but
                        #   we still wait for prior process to finish.
                        if transfer_processes is not None:
                            self.log.info(f"Waiting for {filetype} transfer process "
                                          "to complete.")
                            for p in transfer_processes:
                                p.join()
                        if img_storage_dir is not None:
                            filepath_dests = [img_storage_dir / f for f in filenames]
                            self.log.info("Starting transfer process for "
                                          f"{filepath_dests}.")
                            transfer_processes = [DataTransfer(srcs,dests) for srcs, dests in
                                                  zip(filepath_srcs, filepath_dests)]
                            for p in transfer_processes:
                                p.start()

                        self.tiles_acquired += 1
                        self.tile_time_s = time() - tile_start
                    self.stage_x_pos += x_grid_step_um * STEPS_PER_UM
                self.stage_y_pos += y_grid_step_um * STEPS_PER_UM

        finally:
            if transfer_processes is not None:
                self.log.info("Joining file transfer processes.")
                for p in transfer_processes:
                    p.join()
            if not self.overview_set.is_set():
                dest = str(img_storage_dir) if img_storage_dir != None else str(local_storage_dir)
                for img in self.overview_imgs:
                    DataTransfer(Path(img), Path(dest[:dest.find('\micr')]+img[img.find('\overview_img_'):])).start()

            self.log.info(f"Closing NI tasks")
            self.ni.stop()
            self.log.info(f"Closing camera")
            self.camera.stop()
            for wl, specs in self.cfg.laser_specs.items():
                if str(wl) in self.lasers:
                    self.lasers[str(wl)].disable()

            # Reset values used in scan
            self.active_lasers = None
            self.total_tiles = None
            self.x_y_tiles = None
            self.est_run_time = None
            self.tiles_acquired = 0
            self.tile_time_s = 0

            # Move back to start position
            self.sample_pose.move_absolute(x=self.start_pos['x'], wait=False)
            self.wait_to_stop('x', self.start_pos['x'])  # wait_to_stop uses SAMPLE POSE
            self.sample_pose.move_absolute(y=self.start_pos['y'], wait=False)
            self.wait_to_stop('y', self.start_pos['y'])
            self.sample_pose.move_absolute(z=self.start_pos['z'], wait=False)
            self.wait_to_stop('z', self.start_pos['z'])
            self.log.info(f'Stage moved to {self.sample_pose.get_position()}')

            self.start_pos = None
            acquisition_params = {'session_end_time': datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
                                  'tags': ['schema']}
            self.log.info("acquisition parameters", extra=acquisition_params)

    def _collect_stacked_tiff(self, slow_scan_axis_position: float,
                              tile_count, tile_spacing_um: float,
                              filepath_srcs: list[Path],
                              filetype: str,
                              acquisition_style: str = 'interleaved'):

        self.log.info(f"Configuring stage scan parameters")
        self.log.info(f"Starting scan at Z = {self.stage_z_pos / STEPS_PER_UM / 1000} mm")
        self.sample_pose.setup_finite_tile_scan('z', 'x',
                                                fast_axis_start_position=self.stage_z_pos / STEPS_PER_UM / 1e3,
                                                slow_axis_start_position=slow_scan_axis_position,
                                                slow_axis_stop_position=slow_scan_axis_position,
                                                tile_count=tile_count, tile_interval_um=tile_spacing_um,
                                                line_count=1) if not self.simulated else print('Setting up tile scan')
        # tile_spacing_um = 0.0055 um (property of stage) x ticks
        # Specify fast axis = Tiger x, slow axis = Tiger y,
        frames = (tile_count * len(self.cfg.imaging_wavelengths)) if acquisition_style == 'interleaved' else tile_count
        if self.overview_set.is_set():
            self.stack = [None] * (tile_count)  # Create buffer the size of stacked image

        self.log.info(f"Configuring Camera")
        self.camera.setup_stack_capture(filepath_srcs[0],frames,filetype)   #FIXME: Don't be a filepath list
        self.camera.start()
        self.ni.start()
        self.log.info(f"Starting scan.")
        self.tigerbox.start_scan() if not self.simulated else print('Started')

        prev_frame_count = 0
        self.latest_frame_layer = 0
        while self.ni.counter_task.read() < tile_count:
            if self.simulated:
                tifffile.imwrite(filepath_srcs[0],np.ones((self.cfg.sensor_column_count,
                                                           self.cfg.sensor_row_count)), append=True, bigtiff=True)
            self.frame = self.camera.grab_frame()
            self.latest_frame_layer = self.camera.grab_frame_count()
            if self.camera.grab_frame_count() != prev_frame_count:
                prev_frame_count = self.camera.grab_frame_count()
                self.log.info(f'Total frames: {frames} '
                              f'-> Frames collected: {self.camera.grab_frame_count()}')
            else:
                print('No new frames')
            sleep(self.cfg.get_period_time() + self.cfg.jitter_time_s) if not self.simulated else sleep(.01)

        self.log.info('NI task completed')
        self.log.info('Stopping NI Card')
        self.ni.stop()
        self.__sim_counter_count = 0
        self.latest_frame_layer = 0     # Resetting frame number to 0 for progress bar in UI

        if self.overview_set.is_set():
            self.create_overview() # If doing an overview image, start down sampling and mips
            self.stack = None  # Clear stack buffer

        self.log.info('Stopping camera')
        self.camera.stop()
        self.log.info('Stack complete')

    def _acquisition_livestream_worker(self):

        """Worker yielding the latest frame and frame id during acquisition"""

        while True:
            if self.frame is not None and self.active_lasers is not None:
                if self.cfg.acquisition_style == 'interleaved' and not self.overview_set.is_set():
                    wl = self.active_lasers[self.latest_frame_layer % (len(self.active_lasers)) - 1]
                else:
                    wl = self.active_lasers[0]
                yield self.frame, wl
            else:
                yield # yield so thread can quit
            sleep(.1)

    def overview_scan(self):

        """Quick overview scan function """

        xtiles, ytiles, self.ztiles = self.get_tile_counts(self.cfg.tile_overlap_x_percent,
                                                           self.cfg.tile_overlap_y_percent,
                                                           .8 * 10,
                                                           self.cfg.volume_x_um,
                                                           300,
                                                           self.cfg.volume_z_um)

        self.overview_stack = None  # Clear previous image overview if any
        self.overview_stack = []  # Create empty array size of tiles
        self.overview_set.set()
        # Y volume is always 1 tile
        self.collect_volumetric_image(self.cfg.volume_x_um, 300,
                                      self.cfg.volume_z_um, .8 * 10,
                                      self.cfg.imaging_wavelengths,
                                      (.8 * 10 / 1000 / ((self.cfg.get_period_time()) + self.cfg.jitter_time_s)),
                                      self.cfg.tile_overlap_x_percent, self.cfg.tile_overlap_y_percent,
                                      self.cfg.tile_prefix, 'Trash', self.cfg.local_storage_dir,
                                      acquisition_style='sequential')

        reshaped_array = [None]*len(self.cfg.imaging_wavelengths)
        for wl in self.cfg.imaging_wavelengths:
            index = self.cfg.imaging_wavelengths.index(wl)
            # Split list of all overview images into seperate channels
            split_image_overview= self.overview_stack[index::len(self.cfg.imaging_wavelengths)] if \
                len(self.cfg.imaging_wavelengths) != 1 else self.overview_stack

            # Create empty array size of overview image
            rows = np.shape(split_image_overview[0])[0]
            cols = split_image_overview[0].shape[1]
            overlap = round((self.cfg.tile_overlap_x_percent / 100) * rows)
            overlap_rows = rows - overlap
            reshaped = np.zeros(((xtiles * rows) - (overlap * (xtiles - 1)), ytiles * cols))

            for x in range(0, xtiles):
                for y in range(0, ytiles):
                    cols = split_image_overview[0].shape[1]  # Account for lost frames
                    if x == xtiles - 1:
                        reshaped[x * overlap_rows:(x * overlap_rows) + rows, y * cols:(y + 1) * cols] = \
                        split_image_overview[
                            0]
                    else:
                        reshaped[x * overlap_rows:(x + 1) * overlap_rows, y * cols:(y + 1) * cols] = \
                        split_image_overview[0][
                        0:overlap_rows]
                    del split_image_overview[0]

            reshaped_array[index] = reshaped

        self.overview_imgs.append(fr'{self.cfg.local_storage_dir}\overview_img_{"_".join(map(str, self.cfg.imaging_wavelengths))}'
                         fr'_{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}.tiff')
        pos = self.sample_pose.get_position()
        tifffile.imwrite(self.overview_imgs[-1],
                         reshaped_array,
                         metadata={'position': {'z': pos['z'], 'x': pos['x'], 'y': pos['y']},
                                    'volume':{'z': self.cfg.volume_z_um, 'x': self.cfg.volume_x_um, 'y': 300},
                                    'tile':{'x': xtiles, 'y': ytiles, 'z': self.ztiles}})

        self.overview_set.clear()
        self.overview_process = None
        self.start_pos = None  # Reset start position


        return reshaped_array, xtiles

    def create_overview(self):

        """Create overview image from a stack"""

        self.stack = [i for i in self.stack if i is not None]           # Remove dropped tiles
        downsampled = [x[0::10, 0::10] for x in self.stack]           # Down sample by 10, scikitimage downscale local mean, gpu downsample
        mipstack = [np.max(x, axis=1) for x in downsampled]             # Max projection
        mipstack = np.array(mipstack)

        # Reshape max
        rows = mipstack[0].shape[0]
        cols = mipstack.shape[0]
        reshaped = np.ones((self.ztiles, rows))
        reshaped[0:cols, :] = mipstack
        self.overview_stack.append(np.rot90(np.array(reshaped)))

    def start_livestream(self, wavelength: list, scout_mode: bool):
        """Repeatedly play the daq waveforms and buffer incoming images."""

        self.setting_up_livestream = True
        # Bail early if it's started.
        if self.livestream_enabled.is_set():
            self.log.warning("Not starting. Livestream is already running.")
            return
        self.log.debug("Starting livestream.")
        self.log.warning(f"Turning on the {wavelength}[nm] lasers.")
        self.scout_mode = scout_mode
        self.setup_imaging_for_laser(wavelength, True)
        self.camera.setup_stack_capture(self.cfg.local_storage_dir, 1000000, 'Trash')
        self.livestream_enabled.set()
        # Launch thread for picking up camera images.
        self.setting_up_livestream = False

    def stop_livestream(self, wait: bool = False):

        self.setting_up_livestream = True
        # Bail early if it's already stopped.
        if not self.livestream_enabled.is_set():
            self.log.warning("Not starting. Livestream is already stopped.")
            return
        wait_cond = "" if wait else "not "
        self.log.debug(f"Disabling livestream and {wait_cond}waiting.")
        self.camera.stop()

        self.ni.stop()
        self.ni.close()

        for laser in self.active_lasers: self.lasers[str(laser)].disable()
        self.active_lasers = None
        self.scout_mode = False
        self.livestream_enabled.clear()
        self.setting_up_livestream = False

    def _livestream_worker(self):
        """Pulls images from the camera and puts them into the ring buffer."""

        self.camera.start()
        self.ni.start()
        if self.scout_mode:
            sleep(self.cfg.get_period_time())
            self.ni.stop()
        self.active_lasers.sort()
        while self.livestream_enabled.is_set():
            if self.simulated:
                sleep(1 / 16)
                blank = np.zeros((self.cfg.sensor_row_count,
                                  self.cfg.sensor_column_count),
                                 dtype=self.cfg.image_dtype)
                noise = np.random.normal(0, .1, blank.shape)
                yield noise + blank, 1

            else:
                try:
                    layer_num =self.camera.grab_frame_count() % (len(self.active_lasers)) - 1 if len(self.active_lasers) > 1 \
                        else -1
                    yield self.camera.grab_frame(), self.active_lasers[layer_num + 1]
                except:
                    yield

            sleep((1 / self.cfg.daq_obj_kwds['livestream_frequency_hz']))

    def setup_imaging_for_laser(self, wavelength: list, live: bool = False):
        """Configure system to image with the desired laser wavelength.
        """
        # Bail early if this laser is already setup in the previously set mode.
        if self.active_lasers == wavelength:
            self.log.info("Skipping daq setup. Laser already provisioned.")
            return
        live_status_msg = " in live mode" if self.livestream_enabled.is_set() else ""
        self.log.info(f"Configuring {wavelength}[nm] laser{live_status_msg}.")

        if self.active_lasers is not None:
            for laser in self.active_lasers: self.lasers[str(laser)].disable()

        if self.cfg.acquisition_style == 'sequential':
            fw_index = self.cfg.laser_specs[str(wavelength[0])]['filter_index']  #TODO: This is a hack for getting wavelength
            self.log.info(f"Setting filter wheel to index {fw_index}")
            with self.stage_query_lock:
                self.filter_wheel.set_index(fw_index)

        # Reprovision the DAQ.
        self._setup_waveform_hardware(wavelength, live)
        self.active_lasers = wavelength
        for laser in self.active_lasers: self.lasers[str(laser)].enable()

    def set_scan_start(self, start: dict):

        """Set start position of scan in sample pose.
        :param start: start position of scan in 1/10 um"""

        self.start_pos = start
        self.log.info(f'Scan start position set to {self.start_pos}')

    def calculate_normalized_dct_shannon_entropy(self, image):
        cPSFSupportDiameter = 3
        return normalized_dct_shannon_entropy.compute(image, cPSFSupportDiameter)

    def close(self):
        """Safely close all open hardware connections."""
        self.camera.close()
        self.ni.close()
        for wavelength, laser in self.lasers.items():
            self.log.info(f"Powering down {wavelength}[nm] laser.")
            laser.disable()
        super().close()
