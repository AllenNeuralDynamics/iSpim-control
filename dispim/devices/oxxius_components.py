from oxxius_laser import OxxiusLaser, OXXIUS_COM_SETUP
from serial import Serial

OXXIUS_COM_SETUP = OXXIUS_COM_SETUP

class LaserHub(OxxiusLaser):
    def __init__(self, prefix: str, serial: Serial):
        self.prefix = prefix
        super().__init__(self, serial)
    def _send(self, msg):
        prefix_msg = f'{self.prefix} {msg}'
        return super()._send(prefix_msg)

