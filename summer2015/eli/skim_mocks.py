import astropy.io 
from astropy.io import fits
import numpy as np
import sys
import random
import math
import scipy
import scipy.spatial
from astropy.cosmology import FlatLambdaCDM

import matplotlib.pylab as plt

#Opening Manera mocks
print 'Reading in mocks...'
infilename = sys.argv[1]
vals = np.loadtxt(infilename,unpack=True,skiprows=4)
ra=vals[0]
dec=vals[1]
z=vals[2]
ipoly=vals[3]
wboss=vals[4]
wcp=vals[5]
wred=vals[6]
veto=vals[7]

plt.figure()
plt.plot(ra,dec,'o',markersize=0.2)

norg = (len(ra))
print "Read in %d galaxies." % (norg)

index = veto==1
index *= (wcp+wred-1)==1
index *= (wboss*wcp*wred)==1

ra = ra[index]
dec = dec[index]
z = z[index]

nselect = (len(ra))
print "Writing out %d galaxies." % (nselect)
print "Kept %.2f" % (nselect/float(norg))

plt.figure()
plt.plot(ra,dec,'o',markersize=0.2)

mockdata=np.column_stack((ra,dec,z))
print len(mockdata)
np.savetxt('default.dat',mockdata)

plt.show()


