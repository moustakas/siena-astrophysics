

import astropy.io 
from astropy.io import fits
import numpy as np
import sys
import random
import math
import scipy
import scipy.spatial
from astropy.cosmology import FlatLambdaCDM
'''
 #Opening FITS file (SDSS Data) 'a'
print 'Reading in FITS Data'
infilename1 = sys.argv[1]
hdulist1=fits.open(infilename1)
hdulist1.info()
h1=hdulist1[1]
data2=h1.data
'''
# Opening txt file (Mocks) 'b'
print 'Reading in Text File'
r=np.loadtxt('cmass_dr10_north_randoms_ir4419.v10.1.release.txt')
ra=r[:,0]
dec=r[:,1]
z=r[:,2]
peculiarv=r[:,3]
weight_cboss=r[:,4]
weight_cp=r[:,5]
weight_z=r[:,6]
veto=r[:,7]
'''
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
'''
# Txt Z Cut
zcut=z>0.43
zcutnew=z[zcut]

zcut1=z<0.7
zcutnew1=z[zcut1]
tot=zcut*zcut1
totrand=r[tot]

ra1=totrand[:,0]
dec1=totrand[:,1]
z1=totrand[:,2]
print 'Distances'
 #Comoving Distances
cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
#comdista=cosmo.comoving_distance(data1['Z'])


comdistb=cosmo.comoving_distance(totrand[:,2])
print 'Radians'

##### Converting RA and Dec to Radians #####

# SDSS
#RArada=(data1['PLUG_RA'])*((math.pi)/180)
#Decrada=((math.pi)/2)-((data1['PLUG_DEC'])*((math.pi)/180))

# Mock
RAradb=(totrand[:,0])*((math.pi)/180)
Decradb=((math.pi)/2)-((totrand[:,1])*((math.pi)/180))
print 'Cartesian'

##### Converting to Cartesian Coordinates #####

# SDSS
#xa=comdista*np.sin(Decrada)*np.cos(RArada)
#ya=comdista*np.sin(Decrada)*np.sin(RArada)
#za=comdista*np.cos(Decrada)


#dr10dat=np.column_stack((data1['PLUG_RA'],data1['PLUG_DEC'],data1['Z'],xa,ya,za))
#np.savetxt('adr10dat.txt',dr10dat)
print 'Cartesian 2'
 #Mock
xb=comdistb*np.sin(Decradb)*np.cos(RAradb)
yb=comdistb*np.sin(Decradb)*np.sin(RAradb)
zb=comdistb*np.cos(Decradb)


mockdat=np.column_stack((ra1,dec1,z1,xb,yb,zb))
np.savetxt('amockdat4419.txt',mockdat)














