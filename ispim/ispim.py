"""Abstraction of the ispim Instrument."""

# FIXME: this should live in the napari gui side, and the function
#   we want to launch in a napari thread should exist as a standalone
#   option here.
from datetime import timedelta, datetime
import calendar
import logging
import numpy as np
from pathlib import Path
from time import perf_counter, sleep, time
from mock import NonCallableMock as Mock
from threading import Thread, Event
from ispim.ispim_config import IspimConfig
from ispim.devices.frame_grabber import FrameGrabber
from ispim.devices.ni import WaveformHardware
from ispim.compute_waveforms import generate_waveforms
from ispim.devices.oxxius_components import LaserHub
from serial import Serial
from tigerasi.tiger_controller import TigerController, STEPS_PER_UM
from tigerasi.device_codes import PiezoControlMode, TTLIn0Mode
from tigerasi.sim_tiger_controller import SimTigerController as SimTiger
from spim_core.spim_base import Spim
from spim_core.devices.tiger_components import SamplePose,FilterWheel
from spim_core.processes.data_transfer import DataTransfer
from obis_laser import ObisLS, LSModulationType, AnalogInputImpedanceType
from vortran_laser.stradus import StradusLaser as Vortran
from vortran_laser.stradus import Cmd, Query
import os
from calliphlox import DeviceState
import cv2

LASER_MODEL_TO_OBJ = \
    {
       "Vortran": Vortran,
        "ObisLS": ObisLS
    }

class Ispim(Spim):

    def __init__(self, config_filepath: str,
                 simulated: bool = False):

        self.log = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        # Log setup is handled in the parent class if we pass in a logger.
        super().__init__(config_filepath, simulated=simulated)

        self.cfg = IspimConfig(config_filepath)
        # Instantiate hardware devices
        self.frame_grabber = FrameGrabber() if not self.simulated else \
            Mock(FrameGrabber)
        self.ni = WaveformHardware(**self.cfg.daq_obj_kwds) if not self.simulated else \
            Mock(WaveformHardware)
        self.tigerbox = TigerController(**self.cfg.tiger_obj_kwds) if not \
            self.simulated else SimTiger(**self.cfg.tiger_obj_kwds)
        self.sample_pose = SamplePose(self.tigerbox, **self.cfg.sample_pose_kwds)
        self.filter_wheel = FilterWheel(self.tigerbox,
                                        **self.cfg.filter_wheel_kwds)

        self.lasers = {}  # populated in _setup_lasers.

        # Extra Internal State attributes for the current image capture
        # sequence. These really only need to persist for logging purposes.
        self.total_tiles = 0  # tiles to be captured.
        self.stage_x_pos = None
        self.stage_y_pos = None

        # camera streams filled in with framegrabber.cameras
        self.stream_ids = [item for item in range(0, len(self.frame_grabber.cameras))] if not self.simulated else [0]

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

        self.latest_frame = None
        self.latest_frame_layer = None

        # start position of scan
        self.start_pos = None
        self.im = None

        self.stack = []
        self.image_overview = None
        self.overview_process = None
        self.overview_set = Event()

    def _setup_camera(self):
        """Configure general settings and set camera settings to those specified in config"""

        #TODO: Intialize offset for shape

        self.frame_grabber.setup_cameras((self.cfg.sensor_column_count,
                                          self.cfg.sensor_row_count))

        # Initializing readout direction of camera(s)
        self.frame_grabber.set_scan_direction(0, self.cfg.scan_direction)

        # Initializing line interval of both cameras
        self.frame_grabber.set_line_interval((self.cfg.exposure_time * 1000000) /
                                             self.cfg.sensor_row_count)

        # Initializing exposure time of both cameras
        # TODO: This is assuming that the line_interval is set the same in
        #  both cameras. Should have some fail safe in case not?
        cpx_line_interval = self.frame_grabber.get_line_interval() if not self.simulated else [15, 15]
        self.frame_grabber.set_exposure_time(self.cfg.slit_width_pix *
                                             cpx_line_interval[0])

    def _setup_lasers(self):
        """Setup lasers that will be used for imaging. Warm them up, etc."""

        for wavelength_str, specs in self.cfg.laser_specs.items():
            if specs['kwds']['port'] =='COMxx':
                self.log.warning(f'Skipping setup for laser {wavelength_str} due to no COM port specified')
                continue
            self.log.debug(f"Setting up {wavelength_str} laser.")
            cls = LASER_MODEL_TO_OBJ[specs['driver']]
            # Mock Lasers. We can't wrap them in Mock(cls(**specs['kwds'])
            laser_obj = cls(**specs['kwds']) if not self.simulated else Mock(cls)
            # Put all lasers in analog mode.
            # TODO: can this be abstracted, not case-specific?
            if cls == ObisLS:
                self.log.debug("Applying analog mode settings to ObisLS laser.")
                laser_obj.set_modulation_mode(LSModulationType.ANALOG)
                laser_obj.set_analog_input_impedance(
                    AnalogInputImpedanceType.TWO_THOUSAND_OHM)
            if cls == Vortran:
                laser_obj.set(Cmd.PulseMode, 0)             # Disable digital modulation
                laser_obj.set(Cmd.FiveSecEmissionDelay, 0)  # Disable five sec delay
                laser_obj.set(Cmd.ExternalPowerControl, 1)  # Enable Analog modulation

            self.lasers[int(wavelength_str)] = laser_obj

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
        #self.tigerbox.set_axis_control_mode(**externally_controlled_axes)

        # TODO, this needs to be buried somewhere else
        # TODO, how to store card # mappings, in config?
        self.tigerbox.set_ttl_pin_modes(in0_mode=TTLIn0Mode.MOVE_TO_NEXT_ABS_POSITION,
                                        card_address=31)

    def _setup_waveform_hardware(self, active_wavelength: list, live: bool = False):

        if not self.livestream_enabled.is_set():       # Only configures daq on the initiation of livestream
            self.log.info("Configuring NIDAQ")
            self.ni.configure(self.cfg.get_period_time(), self.cfg.daq_ao_names_to_channels, len(active_wavelength), live)
        self.log.info("Generating waveforms.")
        _, voltages_t = generate_waveforms(self.cfg, active_wavelength)
        self.log.info("Writing waveforms to hardware.")
        self.ni.assign_waveforms(voltages_t)

    # TODO: this should be a base class thing.s
    def check_ext_disk_space(self, dataset_size):
        self.log.warning("Checking disk space not implemented.")

    def wait_to_stop(self, axis: str, desired_position: int):
        """Wait for stage to stop moving. IN SAMPLE POSE"""
        start = time()
        while self.sample_pose.is_moving():
            pos = self.sample_pose.get_position()
            distance = abs(pos[axis.lower()] - desired_position)
            if distance < 1.0 or time()-start > 60:
                self.tigerbox.halt()
                break
            else:
                self.log.info(f"Stage is still moving! {axis} = {pos[axis.lower()]} -> {desired_position}")
                sleep(0.1)

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
                                      self.cfg.local_storage_dir,
                                      self.img_storage_dir,
                                      self.deriv_storage_dir)

    def collect_volumetric_image(self, volume_x_um: float, volume_y_um: float,
                                 volume_z_um: float,
                                 z_step_size_um : float,
                                 channels: list,
                                 scan_speed_mm_s,
                                 tile_overlap_x_percent: float,
                                 tile_overlap_y_percent: float,
                                 tile_prefix: str,
                                 filetype: str,
                                 local_storage_dir: Path = Path("."),
                                 img_storage_dir: Path = None,
                                 deriv_storage_dir: Path = None, ):
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
        # Check if we will violate storage directory limits with our
        #   run. (Check local and external storage.)
        gigabytes_per_image = self.cfg.bytes_per_image / (1024 ** 3)
        dataset_gigabytes = self.total_tiles * gigabytes_per_image
        try:  # Log errors if they occur.
            self.check_ext_disk_space(dataset_gigabytes)
        except AssertionError as e:
            self.log.error(e)
            raise
        # TODO, test network speeds?
        # TODO, check that networked storage is visible?
        # Log relevant info about this imaging run.
        self.schema_log.info(f'session_start_time, {datetime.now().strftime("%Y-%m-%dT%H:%M:%S")}')
        self.schema_log.info(f'local_storage_directoru, {local_storage_dir}')
        self.schema_log.info(f'external_storage_directory, {img_storage_dir}')
        self.schema_log.info(f'specimen_id,{self.cfg.imaging_specs["subject_id"]}')
        self.schema_log.info(f'subject_id,{self.cfg.imaging_specs["subject_id"]}')
        self.schema_log.info(f'instrument_id, iSpim 1')
        self.schema_log.info(f'chamber_immersion_medium, {self.cfg.immersion_medium}')
        self.schema_log.info(f'chamber_immersion_refractive_index, '
                             f'{self.cfg.immersion_medium_refractive_index}')

        self.log.info(f"Total tiles: {self.total_tiles}.")
        self.log.info(f"Total disk space: {dataset_gigabytes:.2f}[GB].")
        run_time_days = self.total_tiles / (self.cfg.tiles_per_second * 3600. * 24.)
        completion_date = datetime.now() + timedelta(days=run_time_days)
        date_str = completion_date.strftime("%d %b, %Y at %H:%M %p")
        weekday = calendar.day_name[completion_date.weekday()]
        self.log.info(f"Time Esimate: {run_time_days:.2f} days. Imaging run "
                      f"should finish after {weekday}, {date_str}.")
        self.log.info(f"Desired dimensions: {volume_x_um:.1f}[um] x "
                      f"{volume_y_um:.1f}[um] x {volume_z_um:.1f}[um]")
        actual_vol_x_um = self.cfg.tile_size_x_um + (xtiles - 1) * x_grid_step_um
        actual_vol_y_um = self.cfg.tile_size_y_um + (ytiles - 1) * y_grid_step_um
        actual_vol_z_um = z_step_size_um * ((ztiles - 1))
        self.log.info(f"Actual dimensions: {actual_vol_x_um:.1f}[um] x "
                      f"{actual_vol_y_um:.1f}[um] x {actual_vol_z_um:.1f}[um]")
        self.log.info(f"X grid step: {x_grid_step_um} [um]")
        self.log.info(f"Y grid step: {y_grid_step_um} [um]")
        self.log.info(f"Z grid step: {z_step_size_um} [um]")
        self.log.info(f"Z grid step: {z_step_size_um} [um]")
        self.log.info(f'xtiles: {xtiles}, ytiles: {ytiles}, ztiles: {ztiles}')
        # TODO, check if stage homing is necessary?

        # Move sample to preset starting position
        if self.start_pos is not None:
            # Start position is in SAMPLE POSE
            self.log.info(f'Moving to starting position at {self.start_pos["x"]}, '
                          f'{self.start_pos["y"]}, '
                          f'{self.start_pos["z"]}')
            self.sample_pose.move_absolute(x=self.start_pos['x'])
            self.wait_to_stop('x', self.start_pos['x'])  # wait_to_stop uses SAMPLE POSE
            self.sample_pose.move_absolute(y=self.start_pos['y'])
            self.wait_to_stop('y', self.start_pos['y'])
            self.sample_pose.move_absolute(z=self.start_pos['z'])
            self.wait_to_stop('z', self.start_pos['z'])
            self.log.info(f'Stage moved to {self.sample_pose.get_position()}')
            # TODO: If reinstate self.sample_pose.zero_in_place() need to set start_pos back to None
        else:
            self.set_scan_start(self.sample_pose.get_position())

        # Set the sample starting location as the origin.
        # self.sample_pose.zero_in_place()

        transfer_processes = None  # Reference to external tiff transfer process.
        # Stage positions in SAMPLE POSE
        self.stage_x_pos, self.stage_y_pos, self.stage_z_pos = (self.start_pos['x'],
                                                                self.start_pos['y'],
                                                                self.start_pos['z'])
        try:
            for j in range(ytiles):
                # move back to x=0 which maps to z=0
                self.stage_x_pos = self.start_pos['x']  # Both in SAMPLE POSE

                # TODO: handle this through sample pose class, which remaps axes
                self.log.info("Setting speed in Y to 1.0 mm/sec")
                self.tigerbox.set_speed(Z=1.0)  # Z maps to Y

                # TODO: handle this through sample pose class, which remaps axes
                self.log.info(f"Moving to Y = {self.stage_y_pos}.")
                self.tigerbox.move_absolute(z=round(self.stage_y_pos))
                self.wait_to_stop('y', self.stage_y_pos)  # wait_to_stop uses SAMPLE POSE

                for i in range(xtiles):
                    # Move to specified X position
                    # TODO: handle this through sample pose class, which remaps axes
                    self.log.debug("Setting speed in X to 1.0 mm/sec")
                    self.tigerbox.set_speed(Y=1.0)  # Y maps to X
                    self.log.debug(f"Moving to X = {round(self.stage_x_pos)}.")
                    self.tigerbox.move_absolute(y=round(self.stage_x_pos))
                    self.wait_to_stop('x', self.stage_x_pos)  # wait_to_stop uses SAMPLE POSE

                    for channel in channels:

                        # TODO: handle this through sample pose class, which remaps axes
                        # Move to specified Z position
                        self.log.debug("Setting speed in Z to 1.0 mm/sec")
                        self.tigerbox.set_speed(X=1.0)  # X maps to Z
                        self.log.debug("Applying extra move to take out backlash.")
                        z_backup_pos = -STEPS_PER_UM * self.cfg.stage_backlash_reset_dist_um
                        self.tigerbox.move_absolute(x=round(z_backup_pos))
                        self.log.info(f"Moving to Z = {self.stage_z_pos}.")
                        self.tigerbox.move_absolute(x=self.stage_z_pos)
                        self.wait_to_stop('z', self.stage_z_pos)  # wait_to_stop uses SAMPLE POSE

                        self.log.info(f"Setting scan speed in Z to {scan_speed_mm_s} mm/sec.")
                        self.tigerbox.set_speed(X=scan_speed_mm_s)
                        self.log.info(f"Actual speed {self.tigerbox.get_speed('x')}mm/sec.")

                        self.log.info(f"Setting up lasers for active channels: {channel}")
                        self.setup_imaging_for_laser([channel])     # Not changing waveform generator so give a one channel list

                        # Setup capture of next Z stack.
                        filetype_suffix = 'tiff' if filetype == 'Tiff' else 'zarr'  # if filetype is trash, it'll be zarr but doesn't matter

                        #channel_string = '_'.join(map(str, self.active_lasers))
                        filenames = [
                            f"{tile_prefix}_X_{i:0>4d}_Y_{j:0>4d}_Z_{0:0>4d}_ch_{channel}.{filetype_suffix}" #add all channel names
                            for
                            camera in self.stream_ids]

                        os.makedirs(local_storage_dir, exist_ok=True)  # Make local directory if not already created
                        filepath_srcs = [local_storage_dir / f for f in filenames]
                        self.log.info(f"Collecting tile stacks at "
                                      f"({self.stage_x_pos / STEPS_PER_UM}, "
                                      f"{self.stage_y_pos / STEPS_PER_UM}) [um] "
                                      f"for channels {channels} and saving to: {filepath_srcs}")
                        # TODO: consider making z step size a fn parameter instead of
                        #   collected strictly from the config.

                        # Logging for JSON schema
                        self.schema_log.info(f'file_name, {filenames}')
                        self.schema_log.info(f'x_voxel_size, {self.cfg.tile_size_x_um} micrometers')  # size of pixels
                        self.schema_log.info(f'y_voxel_size, {self.cfg.tile_size_y_um} micrometers')
                        self.schema_log.info(f'z_voxel_size, {z_step_size_um} micrometers')
                        self.schema_log.info(f'tile_x_position, {self.stage_x_pos * 0.0001} millimeters')
                        self.schema_log.info(f'tile_y_position, {self.stage_y_pos * 0.0001} millimeters')
                        self.schema_log.info(f'tile_z_positione, {self.stage_z_pos * 0.0001} millimeters')
                        self.schema_log.info(f'lightsheet_angle, 45 degrees')

                        self.schema_log.info(f'channel_name, {channel}')
                        self.schema_log.info(f'laser_wavelength, {channel} nanometers')
                        # laser_power = f'{self.lasers[laser].get(Query.LaserPowerSetting)} milliwatts' if int(
                        #     laser) == 561 else \
                        #     f'{self.lasers[laser].get(Query.LaserCurrentSetting)} percent'  # TODO: convert to mW
                        # self.schema_log.info(f'laser_power: {laser_power}')
                        self.schema_log.info(f'filter_wheel_index: {self.cfg.laser_specs[str(channel)]["filter_index"]}')
                        # Every variable in calculate waveforms
                        for key in self.cfg.laser_specs[str(channel)]['etl']:
                            self.schema_log.info(f'daq etl {key}: {self.cfg.laser_specs[str(channel)]["etl"][key]} volts')
                        for key in self.cfg.laser_specs[str(channel)]['galvo']:
                            self.schema_log.info(f'daq galvo {key}: {self.cfg.laser_specs[str(channel)]["galvo"][key]} volts')

                        # Convert to [mm] units for tigerbox.
                        slow_scan_axis_position = self.stage_x_pos / STEPS_PER_UM / 1000.0
                        self._collect_stacked_tiff(slow_scan_axis_position,
                                                   ztiles,
                                                   z_step_size_um,
                                                   filepath_srcs,
                                                   filetype)

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
                            # ↓TODO, use xcopy transfer for speed
                            transfer_processes = [DataTransfer(filepath_srcs[streams],
                                                               filepath_dests[streams]) for streams in self.stream_ids]
                            for p in transfer_processes:
                                p.start()
                        # TODO, set speed of sample Z / tiger X axis to ~1

                    self.stage_x_pos += x_grid_step_um * STEPS_PER_UM
                self.stage_y_pos += y_grid_step_um * STEPS_PER_UM

        finally:
            # TODO, implement sample pose so below can be uncommented
            # self.log.info("Returning to start position.")
            # self.sample_pose.move_absolute(x=0, y=0, z=0, wait=True)
            if transfer_processes is not None:
                self.log.info("Joining file transfer processes.")
                for p in transfer_processes:
                    p.join()
            self.log.info(f"Closing NI tasks")
            self.ni.close() if not self.overview_set.is_set() else self.ni.stop()   # TODO: Should we just stop ni card anyways in case we want to image after?
            self.log.info(f"Closing camera")
            self.frame_grabber.runtime.abort()
            for wl, specs in self.cfg.laser_specs.items():
                self.lasers[wl].disable()
            self.active_lasers = None
            self.schema_log.info(f'Ending time: {datetime.now().strftime("%Y,%m,%d,%H,%M,%S")}')

    def _collect_stacked_tiff(self, slow_scan_axis_position: float,
                              tile_count, tile_spacing_um: float,
                              filepath_srcs: list[Path],
                              filetype: str):

        self.log.info(f"Configuring stage scan parameters")
        self.log.info(f"Starting scan at Z = {self.stage_z_pos / STEPS_PER_UM / 1000} mm")
        self.sample_pose.setup_finite_tile_scan('z', 'x',
                                                fast_axis_start_position=self.stage_z_pos / STEPS_PER_UM / 1e3,
                                                slow_axis_start_position=slow_scan_axis_position,
                                                slow_axis_stop_position=slow_scan_axis_position,
                                                tile_count=tile_count, tile_interval_um=tile_spacing_um,
                                                line_count=1)
        # tile_spacing_um = 0.0055 um (property of stage) x ticks
        # Specify fast axis = Tiger x, slow axis = Tiger y,


        self.log.info(f"Configuring framegrabber")
        self._setup_camera()
        self.frame_grabber.setup_stack_capture(filepath_srcs,
                                               #(tile_count * len(self.cfg.imaging_wavelengths)),
                                               tile_count,
                                               filetype)
        #TODO: Why doesn't this work?
        # while self.frame_grabber.runtime.get_state() != DeviceState.Armed:  # Check if camera is configured
        #     sleep(.05)

        if self.overview_process is not None:   # If doing an overview image, wait till previous tile is done
            if self.overview_process.is_alive():
                self.overview_process.join()
        self.stack = None  # Clear stack buffer
        self.stack = [None]*(tile_count)  # Create buffer the size of stacked image

        self.frame_grabber.start()
        self.ni.start()
        self.log.info(f"Starting scan.")
        self.tigerbox.start_scan()

        prev_frame_count = 0
        curr_frame_count = 0
        while self.ni.counter_task.read() < tile_count:
            for streams in self.stream_ids:
                frame_count = self.framedata(streams)
            curr_frame_count += frame_count
            if curr_frame_count != prev_frame_count:
                prev_frame_count = curr_frame_count
                self.log.info(f'Total frames: {tile_count} '
                              f'-> Frames collected: {curr_frame_count}')
            else:
                print('No new frames')
            sleep(0.05)
        self.log.info('NI task completed')
        self.ni.stop()
        sleep(5)
        if self.overview_set.is_set():
            self.overview_process = Thread(target=self.create_overview)
            self.overview_process.start()   # If doing an overview image, start down sampling and mips

        self.log.info('Waiting for camera to finish')
        start = time()
        while self.frame_grabber.runtime.get_state() == DeviceState.Running:  # Check if camera is finished
            sleep(.05)
            if time() - start > 10:
                self.log.info('Task timed out')
                break
        self.log.info('Stopping camera')
        self.frame_grabber.runtime.abort()
        self.log.info('Stopping NI Card')
        self.ni.stop()
        self.log.info('Stack complete')

    def _acquisition_livestream_worker(self):

        """Worker yielding the latest frame and frame id during acquisition"""

        while True:
            if self.latest_frame is not None:
                yield self.latest_frame, self.cfg.imaging_wavelengths[self.latest_frame_layer %
                                                                      (len(self.cfg.imaging_wavelengths)) - 1]

            sleep(1/17)



    def framedata(self, stream):

        if a := self.frame_grabber.runtime.get_available_data(stream):
            packet = a.get_frame_count()
            f = next(a.frames())
            self.latest_frame= f.data().squeeze().copy()
            self.latest_frame_layer = f.metadata().frame_id

            if self.overview_set.is_set():
                for f in a.frames():
                    self.stack[f.metadata().frame_id] = f.data().squeeze().copy()

            f = None  # <-- fails to get the last frames if this is held?
            a = None  # <-- fails to get the last frames if this is held?
            logging.debug(
                f"Frames in packet: {packet}"
            )
            return packet
        return 0

    def quick_scan(self):

        """Quick overview scan function """

        xtiles, ytiles, self.ztiles = self.get_tile_counts(self.cfg.tile_overlap_x_percent,
                                                      self.cfg.tile_overlap_y_percent,
                                                      self.cfg.z_step_size_um * 10,
                                                      self.cfg.volume_x_um,
                                                      self.cfg.volume_y_um,
                                                      self.cfg.volume_z_um)

        self.image_overview = None                                      # Clear previous image overview if any
        self.image_overview = []                 # Create empty array size of tiles
        self.overview_set.set()
        self.collect_volumetric_image(self.cfg.volume_x_um, self.cfg.volume_y_um,
                                      self.cfg.volume_z_um,self.cfg.z_step_size_um * 10,
                                      self.cfg.imaging_wavelengths,
                                      ((self.cfg.z_step_size_um * 10/1000) / (((self.cfg.get_daq_cycle_time()+.005) * len(self.cfg.imaging_wavelengths)) +0.01 )),
                                      self.cfg.tile_overlap_x_percent, self.cfg.tile_overlap_y_percent,
                                      self.cfg.tile_prefix,'Trash', self.cfg.local_storage_dir)
        if self.overview_process != None:
            self.overview_process.join()

        # Create empty array size of overview image
        rows = np.shape(self.image_overview[0])[0]
        cols = self.image_overview[0].shape[1]
        #total_columns = self.image_overview[x].shape[1]
        overlap = round((self.cfg.tile_overlap_x_percent / 100) * rows)
        overlap_rows = rows - overlap
        reshaped = np.zeros(((xtiles*rows)-(overlap*(xtiles-1)), ytiles*cols))

        for x in range(0, xtiles):
            for y in range(0, ytiles):
                cols = self.image_overview[0].shape[1] # Account for lost frames
                if x == xtiles-1:
                    reshaped[x * overlap_rows:(x * overlap_rows)+rows, y * cols:(y + 1) * cols] = self.image_overview[0]
                else:
                    reshaped[x * overlap_rows:(x + 1) * overlap_rows, y * cols:(y + 1) * cols] = self.image_overview[0][0:overlap_rows]
                del self.image_overview[0]

        self.overview_set.clear()
        #TODO: How to now overwrite?
        cv2.imwrite(
            fr'{self.cfg.local_storage_dir}\overview_img_{"_".join(map(str, self.cfg.imaging_wavelengths))}.tiff',
            reshaped)  # Save overview
        # Move back to start position
        self.sample_pose.move_absolute(x=self.start_pos['x'])
        self.wait_to_stop('x', self.start_pos['x'])  # wait_to_stop uses SAMPLE POSE
        self.sample_pose.move_absolute(y=self.start_pos['y'])
        self.wait_to_stop('y', self.start_pos['y'])
        self.sample_pose.move_absolute(z=self.start_pos['z'])
        self.wait_to_stop('z', self.start_pos['z'])
        self.log.info(f'Stage moved to {self.sample_pose.get_position()}')

        self.start_pos = None   # Reset start position

        return reshaped, xtiles

    def create_overview(self):

        """Create overview image from a stack"""

        self.stack = [i for i in self.stack if i is not None]           # Remove dropped tiles
        self.stack = np.array(self.stack, dtype=object)
        half_col = round(self.cfg.sensor_column_count/2)
        # only use slitwidth pixels
        slit_width = [x[:,half_col-self.cfg.slit_width_pix:half_col+self.cfg.slit_width_pix] for x in self.stack[0:-1]]
        #downsampled = [x[0::10, 0::10] for x in slit_width]           # Down sample by 10, scikitimage downscale local mean, gpu downsample
        mipstack = [np.max(x, axis=1) for x in slit_width]             # Max projection
        mipstack = np.array(mipstack)

        # Reshape max
        rows = mipstack[0].shape[0]
        cols = mipstack.shape[0]
        reshaped = np.ones((self.ztiles, rows))
        reshaped[0:cols, :] = mipstack
        self.image_overview.append(np.rot90(np.array(reshaped)))


    def start_livestream(self, wavelength: list):
        """Repeatedly play the daq waveforms and buffer incoming images."""

        # Bail early if it's started.
        if self.livestream_enabled.is_set():
            self.log.warning("Not starting. Livestream is already running.")
            return
        self.log.debug("Starting livestream.")
        self.log.warning(f"Turning on the {wavelength}[nm] lasers.")
        self.setup_imaging_for_laser(wavelength, live=True)
        self.frame_grabber.setup_stack_capture([self.cfg.local_storage_dir], 1000000, 'Trash')
        self.livestream_enabled.set()
        # Launch thread for picking up camera images.

    def stop_livestream(self, wait: bool = False):
        # Bail early if it's already stopped.

        if not self.livestream_enabled.is_set():
            self.log.warning("Not starting. Livestream is already stopped.")
            return
        wait_cond = "" if wait else "not "
        self.log.debug(f"Disabling livestream and {wait_cond}waiting.")
        self.livestream_enabled.clear()
        self.frame_grabber.runtime.abort()      # Abort for livestream because total frames are never being met

        self.ni.stop()
        self.ni.close()

        for laser in self.active_lasers: self.lasers[laser].disable()
        self.active_lasers = None

    def _livestream_worker(self):
        """Pulls images from the camera and puts them into the ring buffer."""

        self.frame_grabber.start()
        self.ni.start()
        self.active_lasers.sort()

        while self.livestream_enabled.is_set():
            if self.simulated:
                sleep(1 / 16)
                blank = np.zeros((self.cfg.sensor_row_count,
                                  self.cfg.sensor_column_count),
                                 dtype=self.cfg.image_dtype)
                noise = np.random.normal(0, .1, blank.shape)
                yield noise + blank, 1

            elif packet := self.frame_grabber.runtime.get_available_data(self.stream_ids[0]):
                f = next(packet.frames())
                metadata = f.metadata()

                # TODO: Why does this work?
                layer_num = metadata.frame_id % (len(self.active_lasers)) - 1 if len(self.active_lasers) > 1 else -1
                im = f.data().squeeze().copy()
                f = None
                packet = None
                sleep((1 / self.cfg.daq_obj_kwds[
                    'livestream_frequency_hz']) * .1)

                yield im, self.active_lasers[layer_num + 1]

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
            for laser in self.active_lasers: self.lasers[laser].disable()

        fw_index = self.cfg.laser_specs[str(wavelength[0])]['filter_index']  #TODO: This is a hack for
        self.filter_wheel.set_index(fw_index)

        # Reprovision the DAQ.
        self._setup_waveform_hardware(wavelength, live)
        self.active_lasers = wavelength
        for laser in self.active_lasers: self.lasers[laser].enable()

    def set_scan_start(self, start):

        """Set start position of scan in sample pose.
        :param start: start position of scan"""

        self.start_pos = start
        self.log.info(f'Scan start position set to {self.start_pos}')

    def close(self):
        """Safely close all open hardware connections."""
        self.tigerbox.ser.close()
        self.frame_grabber.close()
        self.ni.close()
        for wavelength, laser in self.lasers.items():
            self.log.info(f"Powering down {wavelength}[nm] laser.")
            laser.disable()
        super().close()
