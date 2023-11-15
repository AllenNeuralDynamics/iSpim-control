from abc import ABC, abstractmethod


class Laser(ABC):

    @abstractmethod
    def set_setpoint(self):
        pass

    @abstractmethod
    def get_setpoint(self):
        pass

    @abstractmethod
    def get_max_setpoint(self):
        pass

    @abstractmethod
    def disable_cdrh(self):
        pass

    @abstractmethod
    def enable(self):
        pass

    @abstractmethod
    def disable(self):
        pass
