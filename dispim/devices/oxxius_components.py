from oxxius_laser import OxxiusLaser
from serial import Serial

class LaserHub(OxxiusLaser):
    def __init__(self, prefix: str, serial: Serial):
        self.prefix = prefix
        super().__init__(self, serial)

    def _send(self, msg):
        prefix_msg = f'{self.prefix} {msg}'
        return super()._send(prefix_msg)
