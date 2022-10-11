from oxxius_laser import OxxiusLaser
from serial import Serial

class LaserHub(OxxiusLaser):
    def __init__(self, prefix: str, serial: Serial):
        super().__init__(self, serial)
        self.prefix = prefix

    def _send(self, msg):
        prefix_msg = f'{self.prefix} {msg}'
        return super()._send(prefix_msg)

