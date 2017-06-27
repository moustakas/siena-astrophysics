#!/usr/bin/env python
from __future__ import division
import time
import RPi.GPIO as GPIO
import plotly
import datetime
import os
import sys
import json
import requests

LOGGER = True # if not true, wont log!
FIRST_POST = True
spaces = (" "*50)

#init GPIO
GPIO.setmode(GPIO.BCM)
GPIO.setwarnings(False)

def question(question, poss_ans):
    response = False
    while response == False:
        answer = raw_input(question)
        if answer not in poss_ans:
            print b.ERR+'Not a Valid Answer, please try again...'+b.END
        else:
            response = True
            return answer

#Set Colors for the Prompt
class b:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WARN = '\033[93m'
    ERR = '\033[91m'
    END = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.BLUE = ''
        self.GREEN = ''
        self.WARN = ''
        self.FAIL = ''
        self.END = ''

os.system('clear')
print(b.WARN+"========================================================")
print("Plot.ly + Raspberry Pi + TMP36 Temperature Visualization")
print("========================================================\n"+b.END)

#Questions
FILENAME = raw_input(b.BLUE+"Enter your Graph Name: \n>>  "+b.END)
USERNAME = raw_input(b.BLUE+"Enter your Plotly Username: \n>>  "+b.END)
API_KEY = raw_input(b.BLUE+"Enter your Plotly Api Key: \n>>  "+b.END)
STREAM_TOKEN = raw_input(b.BLUE+"Enter your Plotly Stream Token: \n>>  "+b.END)
TEMP_TYPE = question(b.BLUE+"Do you want to plot Celius for Farenheit? (c/f)  \n>>  "+b.END, ['c','f'])
BAR_SCATTER = question(b.BLUE+"Scatter or Bar Chart? (s/b)  \n>>  "+b.END, ['s','b'])
DELAY = int(raw_input(b.BLUE+"How frequently do you want to post data? (in seconds, minimum is 1)  \n>>  "+b.END))
if DELAY < 1:
    print b.ERR+'You chose a frequency less than 1 seconds.... defaulting to 1 seconds'+b.END
    DELAY = 1


#Set Temp Scale for legend
if TEMP_TYPE == 'c':
    TEMP_LEGEND = 'Temperature (C)'
elif TEMP_TYPE == 'f':
    TEMP_LEGEND = 'Temperature (F)'

if BAR_SCATTER == 's':
    chart_type = 'scatter'
elif BAR_SCATTER == 'b':
    chart_type = 'bar'

#Set Layout
LAYOUT = {
    'title' : FILENAME,
    'showlegend' : True,
    'yaxis' : {
        'title' : TEMP_LEGEND
    }
}

# Function to read data from Analog Pin 0 from MCP3008 (don't need to edit)
def readadc(adcnum, clockpin, mosipin, misopin, cspin):
        if ((adcnum > 7) or (adcnum < 0)):
                return -1
        GPIO.output(cspin, True)

        GPIO.output(clockpin, False)  # start clock low
        GPIO.output(cspin, False)     # bring CS low

        commandout = adcnum
        commandout |= 0x18  # start bit + single-ended bit
        commandout <<= 3    # we only need to send 5 bits here
        for i in range(5):
                if (commandout & 0x80):
                        GPIO.output(mosipin, True)
                else:
                        GPIO.output(mosipin, False)
                commandout <<= 1
                GPIO.output(clockpin, True)
                GPIO.output(clockpin, False)

        adcout = 0
        # read in one empty bit, one null bit and 10 ADC bits
        for i in range(12):
                GPIO.output(clockpin, True)
                GPIO.output(clockpin, False)
                adcout <<= 1
                if (GPIO.input(misopin)):
                        adcout |= 0x1

        GPIO.output(cspin, True)

        adcout /= 2       # first bit is 'null' so drop it
        return adcout

# change these as desired - they're the pins connected from the
# SPI port on the ADC to the GPIO Pins on the Raspi
SPICLK = 18
SPIMISO = 23
SPIMOSI = 24
SPICS = 25

# set up the SPI interface pins
GPIO.setup(SPIMOSI, GPIO.OUT)
GPIO.setup(SPIMISO, GPIO.IN)
GPIO.setup(SPICLK, GPIO.OUT)
GPIO.setup(SPICS, GPIO.OUT)
# temperature sensor connected channel 0 of mcp3008
adcnum = 0

#init plotly
py = plotly.plotly(USERNAME, API_KEY)

response = py.plot([{'x': [], 'y': [], 'type': chart_type, 'name': TEMP_LEGEND, 'stream': {'token': STREAM_TOKEN, 'maxpoints': 1000}}], filename=FILENAME, fileopt='extend', layout=LAYOUT)
print response[u'url']

os.system('clear')
print(b.GREEN+"========================================================")
print("Successfully Posted to Plot.ly! Here is your Graph info:")
print("========================================================\n")
print(b.WARN+"Graph URL:       "+response[u'url'])
print("Graph Title:     "+response[u'filename']+'\n'+b.END)
print(b.GREEN+"========================================================\n"+b.END)

def movingAverage(newx, xarray, winlength):
    xarray.append(newx)
    if len(xarray) > winlength:
        xarray.pop(0)
    return sum(xarray) / len(xarray)

def streamdata():
    TEMP_ARRAY = []
    WIN_LENGTH = 10 # Set to 1 to remove averaging or another value to optimize smoothing for use case
    while True:
        # read the analog pin (temperature sensor LM36)
        read_adc0 = readadc(adcnum, SPICLK, SPIMOSI, SPIMISO, SPICS)

        # convert analog reading to millivolts = ADC * ( 3300 / 1024 )
        millivolts = read_adc0 * ( 3300.0 / 1024.0)

        # 10 mv per degree
        temp_C = ((millivolts - 100.0) / 10.0) - 40.0

        # convert celsius to fahrenheit
        temp_F = ( temp_C * 9.0 / 5.0 ) + 32

        # remove decimal point from millivolts
        millivolts = "%d" % millivolts

        if TEMP_TYPE == 'c':
            TEMP = temp_C
        elif TEMP_TYPE == 'f':
            TEMP = temp_F

        # Take moving average to apply smoothing
        AVG_T = movingAverage(TEMP, TEMP_ARRAY, WIN_LENGTH)

        # show only one decimal place for temperature
        TEMP_READING = "%.1f" % AVG_T

        date_stamp = datetime.datetime.now()
        data = {
        'x': date_stamp.strftime('%Y-%m-%d %H:%M:%S.%f'),
        'y': TEMP_READING,
        }
        yield json.dumps(data)+'\n'
        sys.stdout.write('Temperature Reading:  '+ TEMP_READING)
        sys.stdout.flush()
        sys.stdout.write("\b" * 30)
        time.sleep(DELAY)


requests.post('http://stream.plot.ly',
    data=streamdata(),
    headers = {'plotly-streamtoken': STREAM_TOKEN, 'connection': 'keep-alive'})
