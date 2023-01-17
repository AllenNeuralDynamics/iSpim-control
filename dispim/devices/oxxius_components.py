from oxxius_laser import OxxiusLaser
from serial import Serial

class LaserHub(OxxiusLaser):
    def __init__(self, serial: Serial, prefix=None):
        self.prefix = prefix
        super().__init__(self, serial)
    def _send(self, msg):
        prefix_msg = f'{self.prefix} {msg}' if self.prefix != None else f'{msg}'
        return super()._send(prefix_msg)

