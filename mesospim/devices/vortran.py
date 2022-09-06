"""Vortran laser class"""

import serial
import logging


class Vortran:

    def __init__(self, port):
        self.log = logging.getLogger(__name__ + "." + self.__class__.__name__)
        self.log.warning("No setup to do on Vortran laser.")

    def enable(self):
        pass

    def disable(self):
        pass