import numpy as np
import math
import scipy.spatial
from astropy.cosmology import FlatLambdaCDM
import astropy.io
import matplotlib.pylab as plt
r=np.loadtxt('cmass_dr10_north_randoms_ir4500.v10.1.release.txt')
ra=r[:,0]
dec=r[:,1]
z=r[:,2]


##### Z Cut #####

zcut=z>0.43
zcutnew=z[zcut]

zcut1=z<0.7
zcutnew1=z[zcut1]
tot=zcut*zcut1
totrand=r[tot]
del tot
del r
'''
##### Weight Cut ######

print 'Cutting Weights'
wc=weight_cboss>0
wcnew=weight_cboss[wc]

wp=weight_cp>0
wpnew=weight_cp[wp]

wz=weight_z>0
wznew=weight_z[wz]

##### Veto Cut #######

print 'Veto Cut'
vc=veto==1
vcnew=veto[vc]

##### Final Cut #######

print 'Finalizing Random Cut'
tot=wc*wp*wz*vc
totrand=r[tot]
'''
########### Distances ################

ngals_for_calculation = 100000

np.random.seed(1)

a=np.arange(0,len(totrand))
np.random.shuffle(a)
sample=totrand[a[0:ngals_for_calculation]]
del a
del totrand
###### Distances (Para and Perp) ########

print 'Conversions'
# Comoving Distances
cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
comdist=cosmo.comoving_distance(sample[:,2])

# Converting RA and Dec to Radians

RArad=(sample[:,0])*((math.pi)/180)
Decrad=((math.pi)/2)-((sample[:,1])*((math.pi)/180))

# Converting to Cartesian Coordinates

x=comdist*np.sin(Decrad)*np.cos(RArad)
y=comdist*np.sin(Decrad)*np.sin(RArad)
z=comdist*np.cos(Decrad)
coordsa=np.column_stack((x,y,z))
del RArad
del Decrad
del x
del y
del z
del comdist
del sample
# Distances

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
ngals = len(coordsa)
chunk_size = 50
nchunks = ngals_for_calculation/chunk_size
nbins=200
rangeval=300


tot_freq = np.zeros((nbins,nbins)) 

for j in xrange(nchunks):
    lo = j*chunk_size
    hi = (j+1)*chunk_size
    print "Performing calculations for %d chunk: %d-%d" % (j,lo,hi)

    paras = []
    perps = []

    for i in range(lo,hi):
        if i!=ngals-1:
            # First compute R_LOS and dR
            print 'r1'
            r1=coordsa[i]
            print 'RLOS'
            R_LOS1 = (r1 + coordsa[i+1:])/2.
            print 'dr'
            dR1 = coordsa[i+1:] - r1
            print 'RLOS mag'
            R_LOS_mag1 = mag(R_LOS1)

            # Dot product
            print 'R Para'
            R_para1 = (dR1[:,0]*R_LOS1[:,0] + dR1[:,1]*R_LOS1[:,1] + dR1[:,2]*R_LOS1[:,2])/R_LOS_mag1
            print 'dR mag'
            dR_mag1 = mag(dR1)
            # Make use of the Pythagorean theorem
            print 'R perp'
            R_perp1 = np.sqrt(dR_mag1*dR_mag1 - R_para1*R_para1)
           
            print 'Paras'
            paras += R_para1.tolist()
            print 'Perps'
            perps += R_perp1.tolist()
            
            if i%(chunk_size/4)==0:
                print i

    print 'Histogram'
    hist=plt.hist2d(perps,paras,bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
    print 'Total Frequency'
    tot_freq += hist[0]
    print 'Mirroring'
    # Mirror the negative perps
    hist=plt.hist2d(-1*np.array(perps),paras,bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
    tot_freq += hist[0]


    #print type(hist1[0])
    #frequ1=hist1[0]
    #plt.close()

    print tot_freq
#tot_freq[(nbins/2),(nbins/2)]=0
print 'Final Plot'    
extent = [-rangeval,rangeval, -rangeval,rangeval]
fig = plt.figure()
axes = fig.add_subplot(1,1,1)
print 'Imshow'
ret = axes.imshow(tot_freq,extent=extent,interpolation='nearest') #,origin=origin,cmap=cmap,axes=axes,aspect=aspect
plt.show()
np.savetxt('RRtest2d3.txt',tot_freq)
'''
ngals = len(coords)

paras = []
perps = []
nperps = []
for i,r1 in enumerate(coords):
    if i!=ngals-1:
In [16]: run eli/twopoint.py
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
hist=plt.hist2d(newperps,newparas,bins=200,range=((-150,0),(0,150)))
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
#np.savetxt('RRtestseg1.txt',frequ)
'''
