"""Class for NI Hardware."""
import logging


class WaveformHardware:

    def __init__(self):
        self.log = logging.getLogger(__name__ + "." + self.__class__.__name__)

    def configure(self, period_time: float, ao_names_to_channels: dict):
        """Create internal ni tasks and apply them to hardware."""
        pass

    def assign_waveforms(self, voltages_t):
        """write output waveforms to hardware.

        """
        # Waveform order is the order in which the task channels were created.
        # This is just the ao_names_to_channels iteration order.
        pass

    def start(self):
        pass

    def wait_until_done(self):
        pass

    def stop(self):
        pass

    def close(self):
        pass