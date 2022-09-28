from argparse import ArgumentParser
from dispim.dispim_config import DispimConfig
from tigerasi.tiger_controller import TigerController, UM_TO_STEPS
from tigerasi.sim_tiger_controller import TigerController as SimTiger
# TODO: consolidate these later.
from mesospim.spim_base import Spim
from mesospim.devices.tiger_components import SamplePose

if __name__ == "__main__":
    # Argparse for a config file.
    parser = ArgumentParser()
    parser.add_argument("--config_path", type=str, default="config.toml")
    parser.add_argument("--active_wavelength", type=int, default=488)
    # grab a config filepath.
    args = parser.parse_args()
    config = DispimConfig(args.config_path)

    tigerbox = TigerController(**config.tiger_obj_kwds)
    sample_pose = SamplePose(tigerbox, **config.sample_pose_kwds)

    tigerbox.scanr(X=0, Y=1, Z=32)
    tigerbox.scanv(X=0, Y=0 ,Z=1)

    tigerbox.scan()
    # sample_pose.move_absolute(z=1000, wait=True)
    # tigerbox.move_axes_absolute(x=100, wait_for_output=True, wait_for_reply=True)