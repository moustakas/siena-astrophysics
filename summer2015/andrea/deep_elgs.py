#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from sklearn.naive_bayes import GaussianNB

#reading in OII fitsfile and the values for rz and gr
deep2dir = '/home/desi2/deep2/'
oii = fits.getdata(deep2dir + 'deep2egs-oii.fits.gz',1)
rz_oii = oii['CFHTLS_R']-oii['CFHTLS_Z']
gr_oii = oii['CFHTLS_G']-oii['CFHTLS_R']

#finding the galaxies that have a redshift higher than 1 and a concentration of
#OII larger than stated value
zcut = (oii['z']>1.0)*1
oiicut = (oii['oii_3727']>8E-17)*1

#creating variables for the parameters needed for plt.scatter
colors_oii = np.vstack((rz_oii,gr_oii)).T
labels_oii = zcut & oiicut



stars = fits.getdata('deep2egs-stars.fits.gz',1)
rz_stars = stars['CFHTLS_R']-stars['CFHTLS_Z']
gr_stars = stars['CFHTLS_G']-stars['CFHTLS_R']

gnb = GaussianNB()
dofit = gnb.fit(colors_oii,labels_oii)
hh = 0.02
#setting the min and max values
minoii = 8E-17 # erg/s/cm2
rzmin = -0.5
rzmax = 2.0
grmin = -0.5
grmax = 1.5
xx,yy = np.meshgrid(np.arange(rzmin,rzmax,hh),np.arange(grmin,grmax,hh))
zz = dofit.predict(np.c_[xx.ravel(),yy.ravel()])
zz = zz.reshape(xx.shape)


fig = plt.figure()
ax = fig.add_subplot(111)
ax.contour(xx, yy, zz, [0.5], lw = 2,colors='g')
ax.scatter(colors_oii[:,0],colors_oii[:,1],c=labels_oii,marker='o', alpha=0.6)
plt.axis('tight')
plt.xlim([rzmin,rzmax])
plt.ylim([grmin,grmax])
plt.savefig('/home/desi2/deep2_gauss_ex.png')


# 1. make a plot of g-r vs r-z using cfhtls_g, cfhtls_r, ... for both
# the galaxies and the stars - code the stars, and the strong
# [OII]-emitting galaxies differently

!display /home/desi2/deep2_gauss_ex.png



