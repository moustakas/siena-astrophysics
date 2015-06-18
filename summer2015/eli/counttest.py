import numpy as np
import math
import scipy.spatial
from astropy.cosmology import FlatLambdaCDM
import astropy.io 
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

##### Z Cut #####

zcut=z>0.43
zcutnew=z[zcut]

zcut1=z<0.7
zcutnew1=z[zcut1]
tot=zcut*zcut1
totrand=r[tot]

ngals_for_calculation = 5000

np.random.seed(1)

a=np.arange(0,len(totrand))
np.random.shuffle(a)
sample=totrand[a[0:ngals_for_calculation]]
###### Distances (Para and Perp) ########

print 'Conversions'
# Comoving Distances
cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
comdist=cosmo.comoving_distance(sample[:,2])

RArad=(sample[:,0])*((math.pi)/180)
Decrad=((math.pi)/2)-((sample[:,1])*((math.pi)/180))

# Converting to Cartesian Coordinates

x=comdist*np.sin(Decrad)*np.cos(RArad)
y=comdist*np.sin(Decrad)*np.sin(RArad)
z=comdist*np.cos(Decrad)
coordsa=np.column_stack((x,y,z))

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
nbins=200
rangeval=200
frequencies = []

tot_freq = np.zeros((nbins,nbins)) 

for j in xrange(nchunks):
    lo = j*chunk_size
    hi = (j+1)*chunk_size
    print "Performing calculations for %d chunk: %d-%d" % (j,lo,hi)
    print lo
    paras = []
    perps = []

    for i in range(lo,hi):
        if i!=ngals-1:
            # First compute R_LOS and dR
            r1=coordsa[i]
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
    print len(perps)
    #newperps1=np.concatenate((perps1,nperps1))
    #newparas1=np.concatenate((paras1,paras1))

    #print 'Histogram1'

    import matplotlib.pylab as plt
    hist=plt.hist2d(perps,paras,bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
    tot_freq += hist[0]

    # Mirror the negative perps
    hist=plt.hist2d(-1*np.array(perps),paras,bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
    tot_freq += hist[0]


    #print type(hist1[0])
    #frequ1=hist1[0]
    #plt.close()

    print tot_freq
#tot_freq[(nbins/2),(nbins/2)]=0
extent = [-rangeval,rangeval, -rangeval,rangeval]
fig = plt.figure()
axes = fig.add_subplot(1,1,1)
ret = axes.imshow(tot_freq,extent=extent,interpolation='nearest') #,origin=origin,cmap=cmap,axes=axes,aspect=aspect
plt.show()





