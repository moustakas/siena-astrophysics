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
data2=h1.data

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

#North cut
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

# Randomizing a Sample of SDSS Data

a=np.arange(0,len(data1))
np.random.shuffle(a)
samplea=data1[a[0:17500]]


# Randomizing a sample of Mock Data
b=np.arange(0,len(totrand))
np.random.shuffle(b)
sampleb=totrand[b[0:17500]]


##### Finding Values for Spherical Coordinates ####
'''
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
'''
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
for i,r1 in enumerate(coordsa[0:4375]):
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
for i,r1 in enumerate(coordsa[4375:8750]):
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
#############################################################################
paras3 = []
perps3 = []
nperps3 = []
for i,r1 in enumerate(coordsa[8750:13125]):
        # First compute R_LOS and dR
        R_LOS3 = (r1 + coordsb)/2.
        dR3 = coordsb - r1

        R_LOS_mag3 = mag(R_LOS3)

        # Dot product
        R_para3 = (dR3[:,0]*R_LOS3[:,0] + dR3[:,1]*R_LOS3[:,1] + dR3[:,2]*R_LOS3[:,2])/R_LOS_mag3

        dR_mag3 = mag(dR3)
        # Make use of the Pythagorean theorem
        R_perp3 = np.sqrt(dR_mag3*dR_mag3 - R_para3*R_para3)
        negR_perp3 = -1*np.sqrt(dR_mag3*dR_mag3 - R_para3*R_para3)

        paras3 += R_para3.tolist()
        perps3 += R_perp3.tolist()
        nperps3 +=negR_perp3.tolist()
        print i
print len(paras3)
print len(perps3)
newperps3=np.concatenate((perps3,nperps3))
newparas3=np.concatenate((paras3,paras3))

print 'Histogram3'
import matplotlib.pylab as plt
hist3=plt.hist2d(newperps3,newparas3,bins=200,range=((-150,150),(-150,150)))

frequ3=hist3[0]
#############################################################################
paras4 = []
perps4 = []
nperps4 = []
for i,r1 in enumerate(coordsa[13125:len(coordsa)]):
        # First compute R_LOS and dR
        R_LOS4 = (r1 + coordsb)/2.
        dR4 = coordsb - r1

        R_LOS_mag4 = mag(R_LOS4)

        # Dot product
        R_para4 = (dR4[:,0]*R_LOS4[:,0] + dR4[:,1]*R_LOS4[:,1] + dR4[:,2]*R_LOS4[:,2])/R_LOS_mag4

        dR_mag4 = mag(dR4)
        # Make use of the Pythagorean theorem
        R_perp4 = np.sqrt(dR_mag4*dR_mag4 - R_para4*R_para4)
        negR_perp4 = -1*np.sqrt(dR_mag4*dR_mag4 - R_para4*R_para4)

        paras4 += R_para4.tolist()
        perps4 += R_perp4.tolist()
        nperps4 +=negR_perp4.tolist()
        print i
print len(paras4)
print len(perps4)
newperps4=np.concatenate((perps4,nperps4))
newparas4=np.concatenate((paras4,paras4))

print 'Histogram4'
import matplotlib.pylab as plt
hist4=plt.hist2d(newperps4,newparas4,bins=200,range=((-150,150),(-150,150)))

frequ4=hist4[0]

totfrequ=frequ1+frequ2+frequ3+frequ4
'''
############# NOT IN CHUNKS ################################################### 
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

ngals = len(coordsa)

paras = []
perps = []
nperps= []
for i,r1 in enumerate(coordsa):
    # First compute R_LOS and dR
    R_LOS = (r1 + coordsb)/2.
    dR = coordsb - r1

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
###############################################################################

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
np.savetxt('DRtest2d.txt',frequ)

'''
