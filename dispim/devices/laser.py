"""Base class for laser devices so that we have a common interface."""
import abc


class Laser(abc.ABCMeta):

    def enable(self):
        pass

    def disable(self):
        pass