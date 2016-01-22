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
print 'Reading in file...'
infilename = sys.argv[1]

ngals = 10000

if len(sys.argv)>2:
    ngals = int(sys.argv[2])

vals = np.loadtxt(infilename,unpack=True,skiprows=0)
ra=vals[0]
dec=vals[1]
z=vals[2]

norg = len(ra)

index = np.arange(0,len(ra))
print index
np.random.shuffle(index)
print index

index = index[0:ngals]

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

outname = "%s_n%d.dat" % (infilename.split('.dat')[0],ngals)
np.savetxt(outname,mockdata)

plt.show()


