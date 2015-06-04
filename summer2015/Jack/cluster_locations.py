#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astrometry.libkd.spherematch import match_radec
import os

def cluster_matching():
    jackdir = '/home/desi1/repos/siena-astrophysics/summer2015/Jack/fits/'
    imgs = jackdir+'jpgs/'
    fitsfile = '/home/work/decam/release/dr1/decals-bricks-summary.fits'
    grz_brick = fits.getdata(fitsfile,1)
    ww = np.where(grz_brick['has_image_g'] & grz_brick['has_image_r'] & grz_brick['has_image_z'] & grz_brick['has_catalog'])[0]###& grz_brick['nobs_z_max']>=1)[0]
    redmapper = '/global/work/projects/redmapper/redmapper_isedfit_v5.10_centrals.fits.gz'
    rmap = fits.getdata(redmapper,1)

    fig = plt.figure()
    plt.plot(rmap['ra'],rmap['dec'],'ro')
    plt.plot(grz_brick['ra'][ww],grz_brick['dec'][ww],'bo')
    plt.show()
    fig.savefig('redmapper_clusters.png',clobber=True)
    plt.close(fig)

    ### Spherematching the two sets ###

    m1,m2,d12 = match_radec(grz_brick['ra'][ww][:], grz_brick['dec'][ww][:], rmap['ra'][:], rmap['dec'][:], 0.25/2.0)

    print 'Found', len(m1), 'RA,DEC matches'
    grz_out = grz_brick[ww][:][m1]
    rmap_out = rmap[:][m2] ### All redMaPPer clusters in DECalS footprint ###

    rmap_fits = fits.writeto(jackdir+'rmap.fits',rmap_out,clobber=True)
    grz_fits = fits.writeto(jackdir+'grz.fits',grz_out,clobber=True)
    return rmap,rmap_out

cluster_matching()

### Richness Cutout ###

def rich_cutout():
    nrich = 10
    rich = np.argsort(rmap_out['lambda_chisq'])[::-1][0:nrich]
    stamp = '512'

    for ii in range(nrich):
        print(rmap['ra'][rich[ii]], rmap['dec'][rich[ii]])
        ra = '{:6f}'.format(rmap['ra'][rich[ii]])
        dec = '{:6f}'.format(rmap['dec'][rich[ii]])
        outfile = imgs+'cluster_{:06}'.format(rich[ii])
        url = 'http://imagine.legacysurvey.org/jpeg-cutout-decals-dr1?ra='+ra+'&dec='+dec+'&pixscale=2&size='+stamp
        print(url)
        os.system('wget "'+url+'" -O '+outfile+'.jpg')
        os.system("wget "+url.replace('jpeg','fits')+" -O "+outfile+'.fits')
rich_cutout()
