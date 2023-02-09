from serial import Serial
import time

ser = Serial('COM7',baudrate = 9600)
ser.write(b'?IPA\r')
time.sleep(1)
print(ser.in_waiting)
print(ser.read_until(b'\r\n'))
ser.close()
