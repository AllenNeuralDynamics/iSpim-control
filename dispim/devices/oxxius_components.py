from oxxius_laser import OxxiusLaser
from serial import Serial

class LaserHub(OxxiusLaser):
    def __int__(self, prefix: str, serial: Serial):
        super.__init__(serial)
        self.prefix = prefix

    def send(self, msg):
        super._send(f"{self.prefix} {msg}")

