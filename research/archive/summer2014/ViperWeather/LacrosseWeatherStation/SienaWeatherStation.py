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
        windDc = float(data[3])
        windDir = float(data[4])
        gustDc = float(data[5])


        # convert celsius to fahrenheit
        temperature = ( temp_C * 9.0 / 5.0 ) + 32

        #convert wind deci meters to meters
        wind = windDc/10
        gust = gustDc/10

        date_stamp = datetime.datetime.now()

	return date_stamp,temperature,humidity,rain,wind,windDir,gust


ser = serial.Serial('/dev/ttyACM0',9600)
py.sign_in('physuser','aldyw0r26q')

x1=[]
x2=[]
x3=[]
y1=[]
y2=[]
y3=[]
x4=[]
x5=[]
x6=[]
y4=[]
y5=[]
y6=[]

my_data1 = Scatter(x=x1,y=y1, stream=dict(token='tvfuqv0s6g'), name ='Temp(F)')

my_data2 = Scatter(x=x2,y=y2, stream=dict(token='bjo44dghec'),xaxis='x2',yaxis='y2', name = 'Humidity')
#data=Data([my_data1,my_data2])

my_data3 = Scatter(x=x3,y=y3, stream=dict(token='4tv7be960v'),xaxis='x3',yaxis='y3', name = 'Rain')
#data=Data([my_data1,my_data2,my_data3])

my_data4 = Scatter(x=x4,y=y4, stream=dict(token='onmmq775to'),xaxis='x4',yaxis='y4', name = 'Wind Speed')
#data=Data([my_data1,my_data2,my_data3])

my_data5 = Scatter(x=x5,y=y5, stream=dict(token='g97c3g3475'),xaxis='x5',yaxis='y5', name = 'Wind direction')
#data=Data([my_data1,my_data2,my_data3])

my_data6 = Scatter(x=x6,y=y6, stream=dict(token='6l3o4sttj6'),xaxis='x6',yaxis='y6', name = 'Gust')
data=Data([my_data1,my_data2,my_data3,my_data4,my_data5,my_data6])

layout=Layout(title='Siena Weather',xaxis1=XAxis(anchor='y1',title='time',showline=True),yaxis=YAxis(domain=[0,.12],title='Temperature(Degrees F)',showline=True),xaxis2=XAxis(anchor='y2',title='time',showline=True),yaxis2=YAxis(domain=[0.18,0.30],title='Relative Humidity(%)',showline=True),xaxis3=XAxis(anchor='y3',title='time',showline=True),yaxis3=YAxis(domain=[0.36,0.48],title='Rain',showline=True),xaxis4=XAxis(anchor='y4',title='time',showline=True),yaxis4=YAxis(domain=[0.53,0.65],title='Wind Speed(m/s)',showline=True),xaxis5=XAxis(anchor='y5',title='time',showline=True),yaxis5=YAxis(domain=[0.71,0.83],title='Wind Direction(Degrees)',showline=True),xaxis6=XAxis(anchor='y6',title='time',showline=True),yaxis6=YAxis(domain=[0.89,1],title='Gust(m/s)',showline=True), width=900, height = 1200, autosize=False)


my_fig1=Figure(data=data,layout=layout)
py.plot(my_fig1, auto_open = False)

s1 = py.Stream('tvfuqv0s6g')
s2 = py.Stream('bjo44dghec')
s3 = py.Stream('4tv7be960v')
s4 = py.Stream('onmmq775to')
s5 = py.Stream('g97c3g3475')
s6 = py.Stream('6l3o4sttj6')




s1.open()
s2.open()
s3.open()
s4.open()
s5.open()
s6.open()

while True: #while loop for temperature
	pt = getpoint(ser)
	s1.write(dict(x=pt[0], y=pt[1] ))
	s2.write(dict(x=pt[0], y=pt[2]))
        s3.write(dict(x=pt[0], y=pt[3]))
        s4.write(dict(x=pt[0], y=pt[4] ))
        s5.write(dict(x=pt[0], y=pt[5]))
        s6.write(dict(x=pt[0], y=pt[6]))


        time.sleep(5)


s1.close()
s2.close()
s3.close()
s4.close()
s5.close()
s6.close()


