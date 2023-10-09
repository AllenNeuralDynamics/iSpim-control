from oxxius_laser import OxxiusLaser
from ispim.devices.laser import Laser

class OxxiusDriver(OxxiusLaser,Laser):

    def __init__(self, port,
                 prefix=None,
                 intensity_mode = 'current',
                 modulation_mode: str = None,
                 laser_driver_control_mode: str = None,
                 external_control: str = None):
        super(OxxiusDriver, self).__init__(port, prefix, intensity_mode,
                                           modulation_mode, laser_driver_control_mode,
                                           external_control)

