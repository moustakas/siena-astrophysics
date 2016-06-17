#!/usr/bin/env python

"""
Search for Planet 9 in DECaLS/DR2.

Katie Hoag
2016 June 9
Siena College

"""

import os
import numpy as np
from glob import glob

from astropy.io import fits
from astropy.table import vstack, Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astrometry.libkd.spherematch import match_radec


def get_candidates(cat, gfaint=None):
    
    """ This script will select candidates for Planet 9 from the DECaLS DR2
    Tractor catalogs.
    
    """

    # The sigma checker for g, r, and z filters.
    det_g = (cat['decam_flux'][:, 1]*np.sqrt(cat['decam_flux_ivar'][:, 1]) > 5)
    det_r = (cat['decam_flux'][:, 2]*np.sqrt(cat['decam_flux_ivar'][:, 2]) > 5)
    no_z = (cat['decam_flux'][:, 4]*np.sqrt(cat['decam_flux_ivar'][:, 4]) < 1)

    # Remove WISE data and run sigma checker on WISE 1 and 2.   
    no_w1 = (cat['wise_flux'][:, 0]*np.sqrt(cat['wise_flux_ivar'][:, 0]) < 5)
    no_w2 = (cat['wise_flux'][:, 1]*np.sqrt(cat['wise_flux_ivar'][:, 1]) < 5)

    # Candidates must comply with the following parameters to be considered.
    #pdb.set_trace()
    good = (det_g*det_r*no_z*no_w1*no_w2*\
      (cat['brick_primary']*1)*\
      (cat['decam_nobs'][:, 1] == 1)*\
      (cat['decam_nobs'][:, 2] == 1)*\
      (cat['decam_nobs'][:, 4] >= 1)*\
      (cat['tycho2inblob']*1 == 0)*\
      (cat['out_of_bounds']*1 == 0)*\
      (cat['decam_anymask'][:, 1] == 0)*\
      (cat['decam_anymask'][:, 2] == 0)*\
      (cat['decam_anymask'][:, 4] == 0)*\
      (cat['decam_fracflux'][:, 1] < 0.1)*\
      (cat['decam_fracflux'][:, 2] < 0.1)*\
      (cat['decam_fracflux'][:, 4] < 0.1)*\
      (cat['decam_fracmasked'][:, 1] < 0.1)*\
      (cat['decam_fracmasked'][:, 2] < 0.1)*\
      (cat['decam_fracmasked'][:, 4] < 0.1))*1 
    #good = (det_g*det_r*no_z*no_w1*no_w2)*1

    return np.where(good)[0]

def main():
    
    """ This script selects possible candidates for Planet Nine from
    Tractor DR2 catalogs.
    It rules out objects that are already identified and puts the
    possible Planet Nine candidates into a .fits file.
    
    """
    
    data_dir = os.path.join(os.environ.get('HOME'), 'candidatesp9')
    outfile = os.path.join(data_dir, 'planet9-dr2-candidates.fits')

    catfiles = glob('/global/work/decam/release/dr2/tractor/*/tractor-3*.fits')
    ncat = len(catfiles)
    
    asteroid_path = os.path.join(data_dir, 'asteroids_decals_dr2.fits')
    known_asteroids = fits.getdata(asteroid_path, 1)
    
    gfaint = 30.0
    nout = 0
    
    for ii, thisfile in enumerate(catfiles):
        print('Reading {}'.format(thisfile))
        cat = Table(fits.getdata(thisfile, 1))
        cand = get_candidates(cat, gfaint=gfaint)
        #pdb.set_trace()
        if len(cand) > 0:
            if nout == 0:
                out = cat[cand]
                nout = len(out)
            else:
                out = vstack((out, cat[cand]))
                nout = len(out)
        if nout > 0:
                print('Number of candidates so far: ', nout)
        else:
            print('No candidates yet.')
            
    print('Number of images checked: ', ncat)
    #out.write(outfile, overwrite=True)

    pdb.set_trace()

    if nout > 0:
        # Match candidate catalog (out) against known asteroids

        m1, m2, distance = match_radec(known_asteroids['RA0'], known_asteroids['DEC0'], out['ra'], out['dec'], 1.0/3600.0) # matches within 1 arcsecond (good?)
        keep = np.delete(np.arange(nout), m2)
        finalout = out[keep] #finalout = keep ???
        print("Total Number of Candidates: ", len(finalout))
        print('Writing {}'.format(outfile))
        finalout.write(outfile, overwrite=True)


if __name__ == '__main__':
    main()
