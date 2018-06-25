#!/usr/bin/python

import plotly.plotly as py
import time
import datetime
import serial
#from plotly.graph_objs import Figure, Data,Scatter, Layout

from plotly.graph_objs import *

#reads the data from the sensors and splits the line into the correct variables
def getpoint(ser):
        data = ser.readline().split()
        temp_C = float(data[0])
        humidity = float(data[1])
        rain = float(data[2])
        wind = float(data[3])
        windDir = float(data[4])
        gust = float(data[5])

        # convert celsius to fahrenheit
        temperature = ( temp_C * 9.0 / 5.0 ) + 32
#        temperature = "%.1f" % temperature

        date_stamp = datetime.datetime.now()

	return date_stamp,temperature,humidity,temp_C

# token info
ser = serial.Serial('/dev/ttyACM1',9600)
py.sign_in('physuser','aldyw0r26q')

x1=[]
x2=[]
x3=[]
y1=[]
y2=[]
y3=[]
my_data1 = Scatter(x=x1,y=y1, stream=dict(token='tvfuqv0s6g'))

my_data2 = Scatter(x=x2,y=y2, stream=dict(token='bjo44dghec'),xaxis='x2',yaxis='y2')
data=Data([my_data1,my_data2])

my_data3 = Scatter(x=x3,y=y3, stream=dict(token='4tv7be960v'),xaxis='x3',yaxis='y3')
data=Data([my_data1,my_data2,my_data3])

layout=Layout(yaxis=YAxis(domain=[0,.25]),xaxis2=XAxis(anchor='y2'),yaxis2=YAxis(domain=[0.35,0.60]),xaxis3=XAxis(anchor='y3'),yaxis3=YAxis(domain=[0.70,1]))
#layout=Layout(xaxis1=XAxis(domain=[0.,.45]),xaxis2=XAxis( domain=[.55,1]))


my_fig1=Figure(data=data,layout=layout)
py.plot(my_fig1, auto_open = False)
s1 = py.Stream('tvfuqv0s6g')

#my_data2 = Data([Scatter(x=[],y=[], stream=dict(token='bjo44dghec'))])
#my_fig2=Figure(data=my_data2)
#py.plot(my_fig2, auto_open = False)
s2 = py.Stream('bjo44dghec')
s3=py.Stream('4tv7be960v')



s1.open()
s2.open()
s3.open()
while True: #while loop for temperature
	pt = getpoint(ser)
	s1.write(dict(x=pt[0], y=pt[1] ))
	s2.write(dict(x=pt[0], y=pt[2]))
        s3.write(dict(x=pt[0], y=pt[3]))

        time.sleep(5)


s1.close()
s2.close()
s3.close()

