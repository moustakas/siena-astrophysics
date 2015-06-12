#All Cuts Made
# Argument=twokcmass.fits
import astropy.io 
from astropy.io import fits
import numpy as np
import sys
import random
import math
import scipy
import scipy.spatial
from astropy.cosmology import FlatLambdaCDM

# Opening FITS file (SDSS Data) 'a'
print 'Reading in FITS Data'
infilename1 = sys.argv[1]
hdulist1=fits.open(infilename1)
hdulist1.info()
h1=hdulist1[1]
samplea=h1.data

# Opening txt file (Mocks) 'b'
print 'Reading in Text File'
sampleb=np.loadtxt('twokrand.txt')
ra=sampleb[:,0]
dec=sampleb[:,1]
z=sampleb[:,2]


# Comoving Distances
cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
comdista=cosmo.comoving_distance(samplea['Z'])


comdistb=cosmo.comoving_distance(sampleb[:,2])




##### Converting RA and Dec to Radians #####

# SDSS
RArada=(samplea['PLUG_RA'])*((math.pi)/180)
Decrada=((samplea['PLUG_DEC'])*((math.pi)/180))

# Mock
RAradb=(sampleb[:,0])*((math.pi)/180)
Decradb=((sampleb[:,1])*((math.pi)/180))


##### Converting to Cartesian Coordinates #####

# SDSS
xa=comdista*np.sin(Decrada)*np.cos(RArada)
ya=comdista*np.sin(Decrada)*np.sin(RArada)
za=comdista*np.cos(Decrada)
coordsa=np.column_stack((xa,ya,za))




# Mock
xb=comdistb*np.sin(Decradb)*np.cos(RAradb)
yb=comdistb*np.sin(Decradb)*np.sin(RAradb)
zb=comdistb*np.cos(Decradb)
coordsb=np.column_stack((xb,yb,zb))


print 'Finished with conversions! Now to calculate distances...'





def mag(vec):

    m = None
    # First check if it is an 3xn array of coordinates....
    if type(vec[0])==np.ndarray or type(vec[0])==astropy.units.quantity.Quantity:
        m = np.sqrt(vec[:,0]**2 + vec[:,1]**2 + vec[:,2]**2)
    else:
        # Or if it is just the 3 coordinates.
        m = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

    return m
################################################################################

ngals = len(coordsa)

paras1 = []
perps1 = []
nperps1 = []
for i,r1 in enumerate(coordsa[0:1000]):
        # First compute R_LOS and dR
        R_LOS1 = (r1 + coordsb)/2.
        dR1 = coordsb - r1

        R_LOS_mag1 = mag(R_LOS1)

        # Dot product
        R_para1 = (dR1[:,0]*R_LOS1[:,0] + dR1[:,1]*R_LOS1[:,1] + dR1[:,2]*R_LOS1[:,2])/R_LOS_mag1

        dR_mag1 = mag(dR1)
        # Make use of the Pythagorean theorem
        R_perp1 = np.sqrt(dR_mag1*dR_mag1 - R_para1*R_para1)
        negR_perp1 = -1*np.sqrt(dR_mag1*dR_mag1 - R_para1*R_para1)

        paras1 += R_para1.tolist()
        perps1 += R_perp1.tolist()
        nperps1 +=negR_perp1.tolist()
        print i
print len(paras1)
print len(perps1)
newperps1=np.concatenate((perps1,nperps1))
newparas1=np.concatenate((paras1,paras1))

print 'Histogram1'
import matplotlib.pylab as plt
hist1=plt.hist2d(newperps1,newparas1,bins=200,range=((-150,150),(-150,150)))

frequ1=hist1[0]
#############################################################################
paras2 = []
perps2 = []
nperps2 = []
for i,r1 in enumerate(coordsa[1000:len(coordsa)]):
        # First compute R_LOS and dR
        R_LOS2 = (r1 + coordsb)/2.
        dR2 = coordsb - r1

        R_LOS_mag2 = mag(R_LOS2)

        # Dot product
        R_para2 = (dR2[:,0]*R_LOS2[:,0] + dR2[:,1]*R_LOS2[:,1] + dR2[:,2]*R_LOS2[:,2])/R_LOS_mag2

        dR_mag2 = mag(dR2)
        # Make use of the Pythagorean theorem
        R_perp2 = np.sqrt(dR_mag2*dR_mag2 - R_para2*R_para2)
        negR_perp2 = -1*np.sqrt(dR_mag2*dR_mag2 - R_para2*R_para2)

        paras2 += R_para2.tolist()
        perps2 += R_perp2.tolist()
        nperps2 +=negR_perp2.tolist()
        print i
print len(paras2)
print len(perps2)
newperps2=np.concatenate((perps2,nperps2))
newparas2=np.concatenate((paras2,paras2))

print 'Histogram2'
import matplotlib.pylab as plt
hist2=plt.hist2d(newperps2,newparas2,bins=200,range=((-150,150),(-150,150)))

frequ2=hist2[0]
