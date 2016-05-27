#!/usr/bin/env python
'''
Search for Planet 9 in DECaLS/DR3.

Katie Hoag
2016 May 27
Siena College

'''
import os
import numpy as np
from glob import glob

import pdb

from astropy.io import fits
from astropy.table import vstack
from astropy.coordinates import SkyCoord
from astropy import units as u


def getcandidates(cat, gfaint=None):
    '''Select Planet 9 candidates from a Tractor catalog.'''

    # The sigma checker for g, r, and z filters.  
    det_g = (cat['decam_flux'][:, 1]*np.sqrt(cat['decam_flux_ivar'][:, 1]) > 5)
    det_r = (cat['decam_flux'][:, 2]*np.sqrt(cat['decam_flux_ivar'][:, 2]) > 5)
    no_z = (cat['decam_flux'][:, 4]*np.sqrt(cat['decam_flux_ivar'][:, 4]) < 1)

    # Remove WISE data.  
    # Run sigma checker on WISE 1 and 2.  
    no_w1 = (cat['wise_flux'][:, 0]*np.sqrt(cat['wise_flux_ivar'][:, 0]) < 5)
    no_w2 = (cat['wise_flux'][:, 1]*np.sqrt(cat['wise_flux_ivar'][:, 1]) < 5)

    
    # Candidates must comply with the following parameters to be considered.  

    good = (det_g*det_r*no_z*no_w1*no_w2*
            (cat['brick_primary'] == 'T')*
            (cat['decam_nobs'][:, 1] == 1)*
            (cat['decam_nobs'][:, 2] == 1)*
            (cat['decam_nobs'][:, 4] >= 1)*
            (cat['out_of_bounds'] == 'F')*
            (cat['decam_anymask'][:, 1] == 0)*
            (cat['decam_anymask'][:, 2] == 0)*
            (cat['decam_anymask'][:, 4] == 0)*
            (cat['tycho2inblob'] == 'F')*
            (cat['decam_fracflux'][:, 1] < 0.1)*
            (cat['decam_fracflux'][:, 2] < 0.1)*
            (cat['decam_fracflux'][:, 4] < 0.1)*
            (cat['decam_fracmasked'][:, 1] < 0.1)*
            (cat['decam_fracmasked'][:, 2] < 0.1)*
            (cat['decam_fracmasked'][:, 4] < 0.1))*1
    return np.where(good)[0]
    
    #pdb.set_trace()  # Runs Python Debugger on code up to this line.   

    #wgood = np.where(good)

    # obj_coord = SkyCoord() # known data  
    #cat_coord = SkyCoord(cat['ra'], cat['dec'])  # in degrees  
    #separation = obj_coord.separation(cat_coord)  # Separation might be good to use  
    #if separation <        # small enough for the objects to overlap


def main():

    datadir = os.path.join(os.environ.get('HOME'), 'candidatesp9')
    outfile = os.path.join(datadir, 'planet9-dr3-candidates.fits')

    catfiles = glob(os.path.join(datadir, 'tractor-*.fits'))
    ncat = len(catfiles)

    gfaint = 30.0
    nout = 0
    for ii, thisfile in enumerate(catfiles):
        print('Reading {}'.format(thisfile))
        cat = fits.getdata(thisfile, 1) 
        cand = getcandidates(cat, gfaint=gfaint)

        if len(cand) > 0:
            
            if nout > 0:
                out = cat[cand]
                nout = len(out)
            else:
                out = vstack((out, cat[cand]))
                nout = len(out)

    if nout > 0:
        # Match candidate catalog (out) against known asteroids
        # outcoord = SkyCoord(ra=out['ra'], dec=out['dec'])
        # ...
        # ...
        # finalout = out[non-matching-objects]
        
        print('Writing {}'.format(outfile))
        out.write(outfile, clobber=True)
        #finalout.write(outfile, clobber=True)

if __name__ == '__main__':
    main()
