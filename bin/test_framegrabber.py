import calliphlox
from calliphlox import DeviceKind, Trigger
import time

runtime = calliphlox.Runtime()
dm = runtime.device_manager()
p = runtime.get_configuration()

print(p)
print(p.camera)
print(p.camera.settings)
p.camera.identifier = dm.select(DeviceKind.Camera, name='C15440-20UP')
p.camera.settings.binning = 1
p.camera.settings.shape = (2304, 2304)
p.camera.settings.exposure_time_us = 1725
p.storage.identifier = dm.select(DeviceKind.Storage, name='Tiff')
p.storage.settings.filename = "out.tif"
p.max_frame_count = 10
p.frame_average_count = 0  # disables
runtime.set_configuration(p)
runtime.start()

nframes = 0
while nframes < p.max_frame_count:

    if a := runtime.get_available_data():
        packet = a.get_frame_count()
        f = next(a.frames())
        im = f.data().squeeze()

        f = None  # <-- will fail to get the last frames if this is held?
        a = None  # <-- will fail to get the last frames if this is held?
        nframes += packet
        print(nframes)
        time.sleep(0.1)

runtime.stop()