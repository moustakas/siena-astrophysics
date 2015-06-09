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
print 'Finished reading in data....'

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
sample=data1[a[0:3000]]

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

print 'Finished with conversions! Now to calculate distances...'

# Vectors! (Note that stepcoord is initially the value behind coords[n])  

################################################################################
# Matt's attempt at this.
################################################################################
# Loop over the coordinates and call them r1.
# Then use the np.array functionality for r2 and for 
# calculating all the distances to r1.

################################################################################
# Helper function to calculate the magntiude of a vector.
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

ngals = len(coords)

paras = []
perps = []

for i,r1 in enumerate(coords):
    if i!=ngals-1:
        # First compute R_LOS and dR
        R_LOS = (r1 + coords[i+1:])/2.
        dR = coords[i+1:] - r1

        R_LOS_mag = mag(R_LOS)

        # Dot product
        R_para = (dR[:,0]*R_LOS[:,0] + dR[:,1]*R_LOS[:,1] + dR[:,2]*R_LOS[:,2])/R_LOS_mag

        dR_mag = mag(dR)
        # Make use of the Pythagorean theorem
        R_perp = np.sqrt(dR_mag*dR_mag - R_para*R_para)

        paras += R_para.tolist()
        perps += R_perp.tolist()

        #print R_para,R_perp
        print i
print paras[0:10]
print perps[0:10]
print len(paras)
print len(perps)

################################################################################
print 'Histogram'
import matplotlib.pylab as plt
hist=plt.hist2d(perps,paras,bins=50,range=((0,150),(-150,150)))
#plt.xlabel('Galactic Distances (Mpc)')
#plt.ylabel('Frequency')
#plt.title('Galactic Distance Distribution of 5000 Random CMASS Galaxies')
plt.show()
frequ=hist[0]
distperp=hist[1]
distpar=hist[2]
mperp=distperp[0:-1]
mpar=distpar[0:-1]
diffperp=np.diff(distperp)/2
diffpar=np.diff(distpar)/2
midperp=mperp+diffperp
midpar=mpar+diffpar
vals=np.column_stack((midperp,midpar,frequ))
np.savetxt('DDtest2d.txt',frequ)
