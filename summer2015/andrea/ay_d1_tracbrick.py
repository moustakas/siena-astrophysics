#!/usr/bin/env python


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

dat = fits.getdata('/home/desi2/dr1/decals-bricks.fits',1)
dat.columns

ww = np.where(dat['has_image_g']&dat['has_image_r']&dat['has_image_z']&
dat['has_catalog'])[0]

plt.plot(dat['ra'][ww],dat['dec'][ww],'bo')
plt.savefig('bricks_grz.png')

plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('bricks with the g-r-z')



