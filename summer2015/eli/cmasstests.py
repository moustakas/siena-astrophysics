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
'''
print 'Randomizing A Sample'
a=np.arange(0,len(data1))
np.random.shuffle(a)
sample=data1[a[0:5000]]
'''
ngals_for_calculation = 1000

np.random.seed(1)

a=np.arange(0,len(data1))
np.random.shuffle(a)
sample=data1[a[0:ngals_for_calculation]]
###### Distances (Para and Perp) ########

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
chunk_size = 250
nchunks = ngals_for_calculation/chunk_size

frequencies = []

tot_freq = np.zeros((200,200)) 

for j in xrange(nchunks):
    lo = j*chunk_size
    hi = (j+1)*chunk_size
    print "Performing calculations for %d chunk: %d-%d" % (j,lo,hi)

    paras = []
    perps = []
    
    for i,r1 in enumerate(coordsa[lo:hi]):
        if i!=ngals-1:
            # First compute R_LOS and dR
            R_LOS1 = (r1 + coordsa[i+1:])/2.
            dR1 = coordsa[i+1:] - r1
            R_LOS_mag1 = mag(R_LOS1)

            # Dot product
            R_para1 = (dR1[:,0]*R_LOS1[:,0] + dR1[:,1]*R_LOS1[:,1] + dR1[:,2]*R_LOS1[:,2])/R_LOS_mag1

            dR_mag1 = mag(dR1)
            # Make use of the Pythagorean theorem
            R_perp1 = np.sqrt(dR_mag1*dR_mag1 - R_para1*R_para1)
            #negR_perp1 = -1*R_perp1
            paras += R_para1.tolist()
            perps += R_perp1.tolist()
            #nperps1 += negR_perp1.tolist()
            if i%(chunk_size/4)==0:
                print i

    #print len(paras)
    #print len(perps)
    #newperps1=np.concatenate((perps1,nperps1))
    #newparas1=np.concatenate((paras1,paras1))

    #print 'Histogram1'

    import matplotlib.pylab as plt
    hist=plt.hist2d(perps,paras,bins=200,range=((-150,150),(-150,150)))
    tot_freq += hist[0]

    # Mirror the negative perps
    hist=plt.hist2d(-1*np.array(perps),paras,bins=200,range=((-150,150),(-150,150)))
    tot_freq += hist[0]


    #print type(hist1[0])
    #frequ1=hist1[0]
    #plt.close()

    print tot_freq
tot_freq[100,100]=0
extent = [-150,150, -150,150]
fig = plt.figure()
axes = fig.add_subplot(1,1,1)
ret = axes.imshow(tot_freq,extent=extent,interpolation='nearest') #,origin=origin,cmap=cmap,axes=axes,aspect=aspect
plt.show()
