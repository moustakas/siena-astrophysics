#!/usr/bin/env python

import galsim
import numpy as np
import os
import math
import matplotlib.pyplot as plt
from astropy.io import fits

def decals_stamps():

    nfake = 10
    nstamp = 10
    nx_tiles = 10.0
    ny_tiles = 10.0
    stamp_xsize = 40.0
    stamp_ysize = 40.0
    pixel_scale = 0.2
    random_seed = 123456
    gsparams = galsim.GSParams(maximum_fft_size=2**16)
    
    root_dir = '/home/desi3/decals-stamps/'

    sersicn = np.random.uniform(low=0.8,high=5.0,size=nstamp)
    r50 = np.random.uniform(low=0.1,high=3.0,size=nstamp)
    phi = np.random.uniform(low=0.0,high=180.0,size=nstamp)
    ba = np.random.uniform(low=0.2,high=1.0,size=nstamp)
   
    tbhdu = fits.BinTableHDU.from_columns(
        [fits.Column(name='sersicn',format='f4',array=sersicn),
         fits.Column(name='r50',format='f4',array=r50),
         fits.Column(name='phi',format='f4',array=phi),
         fits.Column(name='ba',format='f4',array=ba)])
    tbhdu.writeto(root_dir+'decals_fake_priors.fits',clobber=True)

    for ii in range(nstamp):
        #print(sersicn[ii],r50[ii],ba[ii],phi[ii])
        gal = galsim.Sersic(n=sersicn[ii],half_light_radius=r50[ii],
                            gsparams=gsparams)
        #psf = galsim.Gaussian(flux=flux, sigma=sigma)
        gal = gal.shear(q=ba[ii],beta=phi[ii]*galsim.degrees)
        image = gal.drawImage()#scale=pixel_scale)
        
        outfile = root_dir+'stamp_{:04}.fits'.format(ii)
        print('Writing '+outfile)
        image.write(outfile,clobber=True)

if __name__ == "__main__":
    decals_stamps()
