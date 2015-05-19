#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astrometry.libkd.spherematch import match_radec

fitsfile = '/home/work/decam/release/dr1/decals-bricks.fits'
grz_brick = fits.getdata(fitsfile,1)
#print dat.columns

ww = np.where(grz_brick['has_image_g'] & grz_brick['has_image_r'] & grz_brick['has_image_z'] & grz_brick['has_catalog'])[0]
#print len(ww)


#plt.plot(dat['ra'][ww],dat['dec'][ww],'bo')
plt.savefig('radec_clusters.png',clobber=True)

redmapper = '/global/work/projects/redmapper/redmapper_isedfit_v5.10_centrals.fits.gz'
rmap = fits.getdata(redmapper,1)
#print dat2.columns

fig = plt.figure()

plt.plot(rmap['ra'],rmap['dec'],'go')
plt.plot(grz_brick['ra'][ww],grz_brick['dec'][ww],'bo')
#this order allows the smaller data set to be in the foreground
plt.show()

fig.savefig('redmapper_clusters.png',clobber=True)


#Spheremapping the two sets

m1,m2,d12 = match_radec(grz_brick['ra'][:], grz_brick['dec'][:], rmap['ra'][:], rmap['dec'][:], 0.25/2.0):

print len(I),'matches'
print 'Found', len(I), 'RA,DEC matches'
col1 = fits.Column(name='Brick_grz', format='E', array=m1)
col2 = fits.Column(name='RedMapper', format='E', array=m2)
cols = fits.ColDefs([col1,col2])
tbhdu = fits.BinTableHDU.from_columns([fits.Column(name='Brick_grz', format='E', array=m1),fits.Column(name='RedMapper', format='E', array=m2)])
merge_fits = tbhdu.writeto('brick_redmapper_merge.fits')

