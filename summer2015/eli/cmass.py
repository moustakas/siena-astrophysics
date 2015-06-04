#### CMASS DATA #####
#### Argument = cmass.fits ####
import astropy.io 
from astropy.io import fits
import numpy as np
import sys
import random
import math
import scipy
from scipy import spatial
infilename1 = sys.argv[1]
hdulist1=fits.open(infilename1)
hdulist1.info()
h1=hdulist1[1]
print 'Reading in Data'
data2=h1.data

#### North Sector Cut ####

Ns1=data2['PLUG_RA']>90
Ns1new=data2['PLUG_RA'][Ns1]
Ns2=data2['PLUG_RA']<280
Ns2new=data2['PLUG_RA'][Ns2]
Ns3=data2['PLUG_DEC']>-10
Ns3new=data2['PLUG_DEC'][Ns3]
Ns4=data2['PLUG_DEC']<80
Ns4new=data2['PLUG_DEC'][Ns4]


tot=Ns1*Ns2*Ns3*Ns4
data1=data2[tot]




###### Distances for Random Galaxies ########
print 'Conversions'
a=np.arange(0,len(data1))
np.random.shuffle(a)
sample=data1[a[0:5000]]
H=70
c=2.99792458*10**5
Zshift=sample['Z'] #Redshift Values
Vr=Zshift*c #Recession Velocity of Galaxy in km/s
rho=Vr/H  #Distance from Earth to Galaxy in Mpc

# Converting RA and Dec to Radians

RArad=(sample['PLUG_RA'])*((math.pi)/180)
Decrad=((sample['PLUG_DEC'])*((math.pi)/180))

# Converting to Cartesian Coordinates

x=rho*np.sin(Decrad)*np.cos(RArad)
y=rho*np.sin(Decrad)*np.cos(RArad)
z=rho*np.cos(Decrad)
coords=np.column_stack((x,y,z))

print 'Distances'
#val=scipy.spatial.distance.cdist(coords,coords)
#vflat = val.flatten()
#d=list(set(vflat))
#d.remove(0.0)
d=scipy.spatial.distance.pdist(coords)

# Histogram
print 'Histogram'
import matplotlib.pylab as plt
hist=plt.hist(d,bins=100,range=(0,20))
plt.xlabel('Galactic Distances (Mpc)')
plt.ylabel('Frequency')
plt.title('Galactic Distance Distribution of 5000 Random CMASS Galaxies')
plt.show()
frequ=hist[0]
dist=hist[1]
m=dist[0:-1]
diff=np.diff(dist)/2
mid=m+diff
vals=np.column_stack((mid,frequ))
np.savetxt('DDtest.txt',vals)
'''
# Scatter

import matplotlib.pylab as plt
plt.plot(data1['PLUG_RA'],data1['PLUG_DEC'],'o',markersize=0.2)
plt.xlabel('Right Ascension (Degrees)')
plt.ylabel('Declination (Degrees)')
plt.title('CMASS Data (DR12)')
#plt.xlim(90,280)
#plt.ylim(-10,80)
plt.show()
'''
