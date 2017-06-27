#!/usr/bin/env python

import os
import numpy as np
#import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import svm
from astroML.utils import completeness_contamination
from astropy.io import fits

filename = '/home/desi2/data'+'/deep2egs-oii.fits.gz'
oii = fits.getdata(filename,1)
rmagcut = (oii['CFHTLS_R']<23.4)*1
oii = oii[:][np.where(rmagcut==1)]
rz = oii['CFHTLS_R']-oii['CFHTLS_Z']
#gr = oii['CFHTLS_G']-oii['CFHTLS_R']
zw1 = oii['CFHTLS_Z']-oii['W1']
#w1w2 = oii['W1']-oii['W2']
zcut = (oii['z']>1.0)*1
oiicut = (oii['oii_3727']>8E-17)*1


filestars = '/home/desi2/data'+'/deep2egs-stars.fits.gz'
stars = fits.getdata(filestars,1)
rz_stars = stars['CFHTLS_R']-stars['CFHTLS_Z']
zw1_stars = stars['CFHTLS_Z']-stars['W1']

#use if you want to classify the rz and zw1 graph
rzw1 = np.vstack((zw1,rz)).T
objtype = zcut & oiicut

rzlim = [-1,2.5]
#grlim = [-0.5,1.5]
w1lim = [-1.5,3.0]
#w2lim = [-0.5,2.0]



    #rz-zw1
plt.figure()
#plt.scatter(rzw1[:,0],rzw1[:,1],c=objtype)
plt.scatter(zw1_stars,rz_stars)
plt.scatter(zw1,rz)

#plt.ylim(rzlim)
#plt.xlim(w1lim)
plt.ylabel('r-z')
plt.xlabel('z-w1')

plt.show()
