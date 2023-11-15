from abc import ABC, abstractmethod


class Camera(ABC):

    @abstractmethod
    def configure(self):
        pass

    @abstractmethod
    def start(self):
        pass

    @abstractmethod
    def stop(self):
        pass

    @abstractmethod
    def grab_frame(self):
        return self.frame

    @abstractmethod
    def grab_frame_count(self):
        return self.frame_count
    @abstractmethod
    def setup_stack_capture(self):
        pass

    @abstractmethod
    def close(self):
        pass

