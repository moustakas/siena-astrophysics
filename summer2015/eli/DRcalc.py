import astropy.io 
from astropy.io import fits
import numpy as np
import sys
import random
import math
import scipy
import scipy.spatial

# Opening FITS file (SDSS Data) 'a'
print 'Reading in FITS Data'
infilename1 = sys.argv[1]
hdulist1=fits.open(infilename1)
hdulist1.info()
h1=hdulist1[1]
data1=h1.data

# Opening txt file (Mocks) 'b'
print 'Reading in Text File'
r=np.loadtxt('cmass_dr10_north_randoms_ir4500.v10.1.release.txt')
ra=r[:,0]
dec=r[:,1]
z=r[:,2]
peculiarv=r[:,3]
weight_cboss=r[:,4]
weight_cp=r[:,5]
weight_z=r[:,6]
veto=r[:,7]
#ztrue=data4560[:,8]
#flag=data4560[:,9]

# Txt Z Cut
zcut=z>0.43
zcutnew=z[zcut]

zcut1=z<0.7
zcutnew1=z[zcut1]
tot=zcut*zcut1
totrand=r[tot]

# Randomizing a Sample of SDSS Data

a=np.arange(0,len(data1))
np.random.shuffle(a)
samplea=data1[a[0:5000]]


# Randomizing a sample of Mock Data
b=np.arange(0,len(totrand))
np.random.shuffle(b)
sampleb=totrand[b[0:20000]]


##### Finding Values for Spherical Coordinates ####

# SDSS
H=70
c=2.99792458*10**5
Zshifta=samplea['Z'] #Redshift Values
Vra=Zshifta*c #Recession Velocity of Galaxy in km/s
rhoa=Vra/H  #Distance from Earth to Galaxy in Mpc

# Mock
Zshiftb=sampleb[:,2] #Redshift Values
Vrb=Zshiftb*c #Recession Velocity of Galaxy in km/s
rhob=Vrb/H  #Distance from Earth to Galaxy in Mpc

##### Converting RA and Dec to Radians #####

# SDSS
RArada=(samplea['PLUG_RA'])*((math.pi)/180)
Decrada=((samplea['PLUG_DEC'])*((math.pi)/180))

# Mock
RAradb=(sampleb[:,0])*((math.pi)/180)
Decradb=((sampleb[:,1])*((math.pi)/180))


##### Converting to Cartesian Coordinates #####

# SDSS
xa=rhoa*np.sin(Decrada)*np.cos(RArada)
ya=rhoa*np.sin(Decrada)*np.cos(RArada)
za=rhoa*np.cos(Decrada)
coordsa=np.column_stack((xa,ya,za))
print coordsa

# Mock
xb=rhob*np.sin(Decradb)*np.cos(RAradb)
yb=rhob*np.sin(Decradb)*np.cos(RAradb)
zb=rhob*np.cos(Decradb)
coordsb=np.column_stack((xb,yb,zb))


##### Distances #####
print 'Distances'
val=scipy.spatial.distance.cdist(coordsa,coordsb)
vflat = val.flatten()
#print "get rid of repeats...."
#d=list(set(vflat)) # Distance Data
#print "Done!"
d = vflat

##### Histogram #####
print 'Histogram'
import matplotlib.pylab as plt
hist=plt.hist(d,bins=100,range=(0,20))
plt.xlabel('Galactic Distances (Mpc)')
plt.ylabel('Frequency')
plt.title('DR with 5000 Galaxies')
plt.show()
frequ=hist[0]
dist=hist[1]
m=dist[0:-1]
diff=np.diff(dist)/2
mid=m+diff
vals=np.column_stack((mid,frequ))
np.savetxt('DRtest.txt',vals)

