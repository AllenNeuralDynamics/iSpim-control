"""Classes for various Mesospim parts, all controlled by an ASI Tigerbox."""

import logging
from time import sleep
from tigerasi.tiger_controller import TigerController


class FilterWheel:
    """Filter Wheel Abstraction from an ASI Tiger Controller."""

    def __init__(self, tigerbox: TigerController, tiger_axis: int = 1):
        """Constructor.

        :param tigerbox: tigerbox hardware.
        :param tiger_axis: which axis the wheel shows up as according to the
            tigerbox.
        """
        self.tigerbox = tigerbox
        self.tiger_axis = tiger_axis #f"FW{wheel_index}"
        self.log = logging.getLogger(__name__ + "." + self.__class__.__name__)

    def get_index(self):
        """return all axes positions as a dict keyed by axis."""
        return self.tigerbox.get_position(str(self.tiger_axis))

    def set_index(self, index: int, wait=True):
        """Set the filterwheel index."""
        cmd_str = f"MP {index}\r\n"
        self.log.debug(f"FW{self.tiger_axis} move to index: {index}.")
        self.tigerbox.send("FW 1\r\n".encode('ascii'))
        self.tigerbox.send(cmd_str.encode('ascii'))
        # TODO: add "busy" check because tigerbox.is_moving() doesn't apply to filter wheels.


class Pose:

    def __init__(self, tigerbox: TigerController):
        """Connect to hardware."""
        self.tigerbox = tigerbox
        self.axes = []
        self.log = logging.getLogger(__name__ + "." + self.__class__.__name__)

    def _move_relative(self, wait: bool = True, **axes: dict):
        axes_moves = "".join([f'{k}={v:.3f} ' for k, v in axes.items()])
        w_text = "" if wait else "NOT "
        self.log.debug(f"Relative move by: {axes_moves} and {w_text}waiting.")
        self.tigerbox.move_axes_relative(**axes, wait_for_output=wait,
                                         wait_for_reply=wait)
        if wait:
            while self.is_moving():
                sleep(0.001)

    def _move_absolute(self, wait: bool = True, **axes: dict):
        """Move the specified axes by their corresponding amounts.

        :param wait: If true, wait for the stage to arrive to the specified
            location. If false, (1) do not wait for the chars to exit the
            serial port, (2) do not wait for stage to respond, and
            (3) do not wait for the stage to arrive at the specified location.
        :param axes: dict, keyed by axis of which axis to move and by how much.
        """
        axes_moves = "".join([f'{k}={v:.3f} ' for k, v in axes.items()])
        w_text = "" if wait else "NOT "
        self.log.debug(f"Absolute move to: {axes_moves}and {w_text}waiting.")
        self.tigerbox.move_axes_absolute(**axes, wait_for_output=wait,
                                         wait_for_reply=wait)
        if wait:
            while self.is_moving():
                sleep(0.001)

    def get_position(self):
        return self.tigerbox.get_position(*self.axes)

    # TODO: do we really need to expose this to the API if we have properly
    #   implemented a wait option?
    def is_moving(self):
        # FIXME: technically, this is true if any tigerbox axis is moving,
        #   but that's all we need for now.
        return self.tigerbox.is_moving()

    def home_in_place(self, *axes):
        """set the specified axes to zero or all as zero if none specified."""
        # We must populate the axes explicitly since the tigerbox is shared
        # between camera stage and sample stage.
        if len(axes) == 0:
            axes = self.axes
        return self.tigerbox.home_in_place(*axes)


class CameraPose(Pose):

    def __init__(self, tigerbox: TigerController):
        super().__init__(tigerbox)
        self.axes = ['N']

    def move_absolute(self, n, wait: bool = True):
        super()._move_absolute(wait, n=n)

    def move_relative(self, n, wait: bool = True):
        axes = {'n': n}
        super()._move_relative(wait, n=n)


class SamplePose(Pose):

    def __init__(self, tigerbox: TigerController):
        super().__init__(tigerbox)
        self.axes = ['X', 'Y', 'Z']

    def move_absolute(self, x=None, y=None, z=None, wait: bool = True):
        # Only specify Non-None axes that we want to move.
        axes = {arg: val for arg, val in locals().items()
                if arg.upper() in self.axes and val is not None}
        super()._move_absolute(wait, **axes)

    def move_relative(self, x=None, y=None, z=None, wait: bool = True):
        # Only specify Non-None axes that we want to move.
        axes = {arg: val for arg, val in locals().items()
                if arg.upper() in self.axes and val is not None}
        super()._move_relative(wait, **axes)
