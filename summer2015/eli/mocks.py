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

a=np.arange(0,len(totrand))
np.random.shuffle(a)
sample=totrand[a[0:10000]]

###### Distances (Para and Perp) ########

print 'Conversions'
# Comoving Distances
cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
comdist=cosmo.comoving_distance(sample[:,2])







'''
H=70
c=2.99792458*10**5
Zshift=sample[:,2] #Redshift Values
Vr=Zshift*c #Recession Velocity of Galaxy in km/s
rho=Vr/H  #Distance from Earth to Galaxy in Mpc
'''
# Converting RA and Dec to Radians

RArad=(sample[:,0])*((math.pi)/180)
Decrad=((sample[:,1])*((math.pi)/180))

# Converting to Cartesian Coordinates

x=comdist*np.sin(Decrad)*np.cos(RArad)
y=comdist*np.sin(Decrad)*np.sin(RArad)
z=comdist*np.cos(Decrad)
coords=np.column_stack((x,y,z))
# Distances
'''
d1=scipy.spatial.distance.pdist(coords)

vflat1 = val1.flatten()
d1=list(set(vflat))
d1.remove(0.0)
'''
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
hist=plt.hist2d(perps,paras,bins=100,range=((0,3000),(-1000,1000)))
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
np.savetxt('RRtest2d.txt',vals)

