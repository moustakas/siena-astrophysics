#! /usr/bin/env python

import os
import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def get_stamp(redshift):
    # Do calculations that depend on redshift
    pixscale = 0.262
    in_dir = os.getenv('HOME')+'/redmapper/' # push to a function
    out_dir = in_dir+'cutouts/'
    rmap = fits.getdata(in_dir+'rmap.fits',1)
    theta = 500.0*cosmo.arcsec_per_kpc_comoving(rmap['z']).value
    stamp = theta/pixscale
    return stamp

def main():
    """Document me
    
    """
    in_dir = os.getenv('HOME')+'/redmapper/' # push to a function
    out_dir = in_dir+'cutouts/'
    
    rmap = fits.getdata(in_dir+'rmap.fits',1)
    
    nrich = 10
    rich = np.argsort(rmap['lambda_chisq'])[::-1][0:nrich]
    print('The top '+str(nrich)+' richest clusters: '+str(rmap['lambda_chisq'][rich]))

    for ii in range(nrich):
        print('Working on cluster {}'.format(ii))
        ra = '{:6f}'.format(rmap['ra'][rich[ii]])
        dec = '{:6f}'.format(rmap['dec'][rich[ii]])
        cluster = 'cluster_{:06}'.format(rich[ii])
        stamp = '{:n}'.format(np.floor(get_stamp(rmap['z'])[rich[ii]])/10*10)

        aa = os.system("$TRACTOR_DIR/projects/desi/runbrick.py --radec "+ra+" "+dec+" --width 600 --height 600 --no-wise --no-sdss --threads 16")
        print(aa)

if __name__ == '__main__':
    main()

