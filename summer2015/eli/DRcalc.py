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
print "Read in text file......"

print "Making cuts......"
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
print "Made cuts....."

# Randomizing a Sample of SDSS Data
ngals_for_calculation = 60000
nrands=75000
np.random.seed(1)

a=np.arange(0,len(data1))
np.random.shuffle(a)
samplea=data1[a[0:ngals_for_calculation]]


# Randomizing a sample of Mock Data
b=np.arange(0,len(totrand))
np.random.shuffle(b)
sampleb=totrand[b[0:nrands]]


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
chunk_size = 1000
nchunks = ngals_for_calculation/chunk_size

frequencies = []

tot_freq = np.zeros((300,300)) 

for j in xrange(nchunks):
    lo = j*chunk_size
    hi = (j+1)*chunk_size
    print "Performing calculations for %d chunk: %d-%d" % (j,lo,hi)

    paras = []
    perps = []
  
    for i,r1 in enumerate(coordsa[lo:hi]):
            # First compute R_LOS and dR
            R_LOS1 = (r1 + coordsb[:])/2
            dR1 = coordsb - r1
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
    #newparas1len=np.concatenate((paras1,paras1))

    #print 'Histogram1'

    import matplotlib.pylab as plt
    hist=plt.hist2d(perps,paras,bins=300,range=((-150,150),(-150,150)))
    tot_freq += hist[0]

    # Mirror the negative perps
    hist=plt.hist2d(-1*np.array(perps),paras,bins=300,range=((-150,150),(-150,150)))
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
np.savetxt('DRtest2d1.txt',tot_freq)

#newperps2=np.concatenate((perps2,nperps2))
#newparas2=np.concatenate((paras2,paras2))

#print 'Histogram2'
#import matplotlib.pylab as plt
#hist2=plt.hist2d(newperps2,newparas2,bins=200,range=((-150,150),(-150,150)))
#
#frequ2=hist2[0]
#plt.close()
#############################################################################
#totfrequ=frequ1+frequ2+frequ3+frequ4
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
