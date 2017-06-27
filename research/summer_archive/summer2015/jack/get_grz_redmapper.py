#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astrometry.libkd.spherematch import match_radec
import os

def main():
    """Document me

    """
    out_dir = os.getenv('HOME')+'/redmapper/'
    redmapper_dir = '/global/work/projects/redmapper/'
    dr1_dir = os.getenv('HOME')+'/dr1/'
    #decals_dir = os.getenv('DECALS_DIR')+'/'

    # read bricks summary file and find bricks with grz
    print('Reading '+dr1_dir+'decals-bricks-summary.fits')
    allbricks = fits.getdata(dr1_dir+'decals-bricks-summary.fits',1)
    these = np.where(allbricks['has_image_g'] & allbricks['has_image_r'] & 
                     allbricks['has_image_z'] & allbricks['has_catalog'] & 
                     (allbricks['nobs_g_med']>2)*1 & (allbricks['nobs_r_med']>2)*1 & 
                     (allbricks['nobs_z_med']>2))[0]
    grz_bricks = allbricks[these]

    print('Reading '+redmapper_dir+'redmapper_isedfit_v5.10_centrals.fits.gz')
    rmap = fits.getdata(redmapper_dir+'redmapper_isedfit_v5.10_centrals.fits.gz',1)

    # Match the two datasets.
    m1, m2, d12 = match_radec(grz_bricks['ra'],grz_bricks['dec'], 
                              rmap['ra'],rmap['dec'],0.25/2.0,nearest=False)

    print('Found {} redMaPPer clusters in DECalS footprint!'.format(len(m1)))
    grz_out = grz_bricks[list(set(m1))]
    rmap_out = rmap[m2]

    rmap_fits = fits.writeto(out_dir+'rmap.fits',rmap_out,clobber=True)
    grz_fits = fits.writeto(out_dir+'grz.fits',grz_out,clobber=True)

    # Make QAplot.  
    # TODO (@thejackparticle): add labels
    fig = plt.figure()
    plt.plot(rmap['ra'],rmap['dec'],'ro')
    plt.plot(grz_bricks['ra'],grz_bricks['dec'],'bo')
    fig.savefig(out_dir+'redmapper_clusters.png',clobber=True)
    plt.close(fig)

if __name__ == '__main__':
    main()
