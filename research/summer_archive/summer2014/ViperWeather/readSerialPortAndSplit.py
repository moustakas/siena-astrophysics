#!/usr/bin/python
import serial
ser = serial.Serial('/dev/ttyACM0',9600)

data = ser.readline().split()
print data
#temp_C = float(data[0])
#humidity = float(data[1])
#rain = float(data[2])
#wind = float(data[3])
#windDir = float(data[4])
#gust = float(data[5])
        
#print temp_C
#,humidity,rain,wind,windDir,gust
