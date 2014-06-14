#!/usr/bin/python

import plotly.plotly as py
import time
import datetime
import serial
from plotly.graph_objs import Figure, Data,Scatter

#reads the data from the sensors and splits the line into the correct variables
def getpoint(ser):
        data = ser.readline().split()
        humidity = float(data[0])
        temp_C = float(data[1])

        # convert celsius to fahrenheit
        temperature = ( temp_C * 9.0 / 5.0 ) + 32
#        temperature = "%.1f" % temperature

        date_stamp = datetime.datetime.now()

	return date_stamp,temperature,humidity

# token info
ser = serial.Serial('/dev/ttyACM1',9600)
py.sign_in('physuser','aldyw0r26q')
my_data1 = Data([Scatter(x=[],y=[], stream=dict(token='tvfuqv0s6g'))])
my_fig1=Figure(data=my_data1)
py.plot(my_fig1, auto_open = False)
s1 = py.Stream('tvfuqv0s6g')

my_data2 = Data([Scatter(x=[],y=[], stream=dict(token='bjo44dghec'))])
my_fig2=Figure(data=my_data2)
py.plot(my_fig2, auto_open = False)
s2 = py.Stream('bjo44dghec')



s1.open()
s2.open()
while True: #while loop for temperature
	pt = getpoint(ser)
	s1.write(dict(x=pt[0], y=pt[1] ))
	s2.write(dict(x=pt[0], y=pt[2]))
        time.sleep(5)


s1.close()
s2.close()

