#! /usr/bin/env python

import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

def get_stamp(redshift):
    # Do calculations that depend on redshift
    pixscale = 0.262
    
    stamp = 512
    return stamp

def main():
    """Document me! (Read PEP8)

    Write me.

    """
    in_dir = os.getenv('HOME')+'/redmapper/' # push to a function
    out_dir = in_dir+'cutouts/'

    rmap = fits.getdata(in_dir+'rmap.fits',1)

    nrich = 10
    rich = np.argsort(rmap['lambda_chisq'])[::-1][0:nrich]
    print(rmap['lambda_chisq'][rich])

    for ii in range(nrich):
        print('Working on cluster {}'.format(ii))
        ra = '{:6f}'.format(rmap['ra'][rich[ii]])
        dec = '{:6f}'.format(rmap['dec'][rich[ii]])
        cluster = 'cluster_{:06}'.format(rich[ii])
        stamp = get_stamp(rmap['z'][rich[ii]])
        
        jpeg_url = 'http://imagine.legacysurvey.org/jpeg-cutout-decals-dr1?ra='+ra+'&dec='+dec+'&pixscale=2&size='+str(stamp)
        fits_url = 'http://imagine.legacysurvey.org/fits-cutout-decals-dr1?ra='+ra+'&dec='+dec+'&pixscale=2&size='+str(stamp)
        print(jpeg_url)
        print(fits_url)
        os.system('wget "'+jpeg_url+'" -O '+out_dir+cluster+'.jpg')
        os.system('wget "'+fits_url+'" -O '+out_dir+cluster+'.fits')

if __name__ == '__main__':
    main()
