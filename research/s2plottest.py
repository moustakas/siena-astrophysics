import s2plot
import numpy as np
from s2plot import *

s2opendo("/s2mono") # opens desired window format

s2swin(-2000,2000,-2000,2000,-2000,2000)  # Axes range (x1,x2,y1,y2,z1,z2)
s2box("BCDET",0,0,"BCDET",0,0,"BCDET",0,0) #Box drawing

#Data
#### CMASS DATA #####
#### Argument = dr10cmassnorth.fits ####
import astropy.io 
from astropy.io import fits
import sys
import random
import math
import scipy
from scipy import spatial
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pylab as plt
import time
infilename1 = sys.argv[1]
hdulist1=fits.open(infilename1)
hdulist1.info()
h1=hdulist1[1]
print 'Reading in Data'
data1=h1.data
print 'Finished reading in data....'

del h1

ngals_for_calculation = 400000
np.random.seed(1)

a=np.arange(0,len(data1))
np.random.shuffle(a)
sample=data1[a[0:ngals_for_calculation]]
###### Distances (Para and Perp) ########
del data1
print 'Conversions'
# Comoving Distances
cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
comdist=cosmo.comoving_distance(sample['Z'])

# Converting RA and Dec to Radians

RArad=(sample['PLUG_RA'])*((math.pi)/180)
Decrad=((math.pi)/2)-((sample['PLUG_DEC'])*((math.pi)/180))

# Coverting to Cartesian Coordinates

x=comdist*np.sin(Decrad)*np.cos(RArad)
y=comdist*np.sin(Decrad)*np.sin(RArad)
z=comdist*np.cos(Decrad)
coordsa=np.column_stack((x,y,z))
del RArad
del Decrad
del sample
del comdist
print max(x)
print min(x)
print max(y)
print min(y)
print min(z)
print max(z)
print len(y)
print len(z)

# Back to S2plot

s2icm("astro",1000,1500)  #color map
s2scir(1000,1500) #indexing 
s2sch(2) # character height
print "labels"
s2lab("X","Y","z","Sample of CMASS data") #labels
print "point draw"
s2pt(400000,x,y,z,1)  #Drawing Points, the first number is the number of points
print "display"
s2disp(-1,1)
print "done"
