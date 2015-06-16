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
ngals_for_calculation = 100000

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
Decrad=((sample['PLUG_DEC'])*((math.pi)/180))

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
chunk_size = 500
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

extent = [-150,150, -150,150]
fig = plt.figure()
axes = fig.add_subplot(1,1,1)
ret = axes.imshow(tot_freq,extent=extent,interpolation='nearest') #,origin=origin,cmap=cmap,axes=axes,aspect=aspect
plt.show()
np.savetxt('DDtest2d1.txt',tot_freq)
'''
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

frequtot=frequ1+frequ2

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
nperps = []
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
        negR_perp = -1*np.sqrt(dR_mag*dR_mag - R_para*R_para)

        paras += R_para.tolist()
        perps += R_perp.tolist()
        nperps +=negR_perp.tolist()
        #print R_para,R_perp
        print i
print paras[0:10]
print perps[0:10]
print len(paras)
print len(perps)
newperps=np.concatenate((perps,nperps))
newparas=np.concatenate((paras,paras))



################################################################################
print 'Histogram'
import matplotlib.pylab as plt
hist=plt.hist2d(newperps,newparas,bins=200,range=((-150,150),(-150,150)))
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
'''
