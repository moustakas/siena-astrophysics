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
from astropy.cosmology import FlatLambdaCDM
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

###### Random Sample #######

print 'Randomizing A Sample'
a=np.arange(0,len(data1))
np.random.shuffle(a)
sample=data1[a[0:200]]

###### Distances (Para and Perp) ########

print 'Conversions'
# Comoving Distances
cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
comdist=cosmo.comoving_distance(sample['Z'])

# Converting RA and Dec to Radians

RArad=(sample['PLUG_RA'])*((math.pi)/180)
Decrad=((sample['PLUG_DEC'])*((math.pi)/180))

# Coverting to Cartesian Coordinates

x=comdist*np.sin(Decrad)*np.cos(RArad)
y=comdist*np.sin(Decrad)*np.sin(RArad)
z=comdist*np.cos(Decrad)
coords=np.column_stack((x,y,z))

# Vectors! (Note that stepcoord is initially the value behind coords[n])  

R_LOS=[]
print 'Line of Sight'
for i in range(0,len(sample)):
    print i
    stepcoord=coords[i]
    for n in range(i+1,len(sample)):
        loscalc=[(stepcoord+coords[n])/2]
        R_LOS.append(loscalc[0])
R_LOS=np.array(R_LOS)

dR=[]
print 'dR'
for i in range(0,len(sample)):
    print i
    stepcoord1=coords[i]
    for n in range(i+1,len(sample)):
        dRcalc=[(coords[n]-stepcoord1)]
        dR.append(dRcalc[0])
dR=np.array(dR)

R_LOSmag=[]
print 'Line of sight magnitude'
for i in range(0,len(R_LOS)):
    print i
    losstep=R_LOS[i]
    losmagcalc=[np.sqrt((losstep[0]**2)+(losstep[1]**2)+(losstep[2]**2))]
    R_LOSmag.append(losmagcalc[0])
R_LOSmag=np.array(R_LOSmag)

rpartop=[]
print 'Numerator of Rpar'
for i in range(0,len(R_LOS)):
    print i
    lstep=R_LOS[i]
    dstep=dR[i]
    rcalc=[(lstep[0]*dstep[0])+(lstep[1]*dstep[1])+(lstep[2]*dstep[2])]
    rpartop.append(rcalc[0])
rpartop=np.array(rpartop)

Rpar=rpartop/R_LOSmag
#Rperp=dR-Rpar   #What?



'''
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
y=rho*np.sin(Decrad)*np.sin(RArad)
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
