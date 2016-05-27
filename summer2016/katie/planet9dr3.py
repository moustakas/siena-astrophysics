#!/usr/bin/env python


# Planet Nine code for DR3

import os
import numpy as np
import glob

datadir = '/home/desi2/candidatesp9/'
known = glob.glob(os.path.join(datadir, 'tractor-*.fits'))

def candidate_bitmask(cat): #brightness????

    # Sigma checker for g, r, and z filters
    det_g = (cat['decam_flux[1, *]']*sqrt(cat['decam_flux_ivar[1, *]']) > 5)
    det_r = (cat['decam_flux[2, *]']*sqrt(cat['decam_flux_ivar[2, *]']) > 5)
    no_z = (cat['decam_flux[4, *]']*sqrt(cat['decam_flux_ivar[4, *]']) < 1)

    # Remove WISE data
    # Sigma checker on WISE 1 and 2
    no_w1 = (cat['wise_flux[0, *]']*sqrt(cat['wise_flux_ivar[0, *]']) < 5)
    no_w2 = (cat['wise_flux[1, *]']*sqrt(cat['wise_flux_ivar[1, *]']) < 5)

    good = np.ndarray(n_elements(cat)) + 1
    
    # Must comply with the following parameters to be considered
    (good and (det_g and det_r and no_z)) * 1
    (good and (cat['brick_primary'] == 'T')) * 1
    (good and (cat['decam_nobs' [1, :]] == 1)) * 1 # g filter
    (good and (cat['decam_nobs' [2, :]] == 1)) * 1 # r filter
    (good and (cat['decam_nobs' [4, :]] > 1)) * 1 # z filter
    (good and (cat['out_of_bounds'] == 'F')) * 1
    (good and (cat['decam_anymask' [1, :]] == 0)) * 1 # g filter
    (good and (cat['decam_anymask' [2, :]] == 0)) * 1 # r filter
    (good and (cat['decam_anymask' [4, :]] == 0)) * 1 # z filter
    (good and (cat['tycho2inblob'] == 'F')) * 1
    (good and (cat['decam_fracflux'[1, :]] < 0.1)) * 1 # g filter
    (good and (cat['decam_fracflux'[2, :]] < 0.1)) * 1 # r filter
    (good and (cat['decam_fracflux'[4, :]] < 0.1)) * 1 # z filter
    (good and (cat['decam_fracmasked'[1, :]] < 0.1)) * 1 # g filter
    (good and (cat['decam_fracmasked'[2, :]] < 0.1)) * 1 # r filter
    (good and (cat['decam_fracmasked'[4, :]] < 0.1)) * 1 # z filter

    good and (no_w1 and no_w2) * 1


    #if # keyword_set(bright)?



   #if np.sum(good) == 0:
    #    return good

    rejects=[]
    wgood = where(good)
    #if cat[wgood].ra == known.ra0 and cat[wgood].dec == known.dec0:
        # spherematch line?
     #   rejects.append(wgood)
      #  return'''
       
    return
