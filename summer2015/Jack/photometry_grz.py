#! /usr/bin/env python

from astropy.table import vstack, Table
from astropy.io import fits
from astrometry.libkd.spherematch import match_radec
import os

def main():
    """Obtain photometry for grz.fits bricks using the tractor fits files which also belong to the redmapper catalog.

    """
    out_dir = os.getenv('HOME')+'/redmapper/'
    in_dir = os.getenv('HOME')+'/dr1/tractor/'
    redmapper_dir = '/global/work/projects/redmapper/'
    decals_dir = os.getenv('DECALS_DIR')+'/'

    grz_fits = fits.getdata(out_dir+'grz.fits',1)
    #grz_fits = grz_fits[:100]
    rmap = fits.getdata(out_dir+'rmap.fits',1)

    cat = None
    for brick in grz_fits:
        brickname = str(brick['brickname'])
        raslice = brickname.replace('tractor-','')[:3]
        print(raslice,brickname)
        catfile = in_dir+raslice+'/tractor-'+brickname+'.fits'
        tractor = fits.getdata(catfile,1)

        m1, m2, d12 = match_radec(tractor['ra'],tractor['dec'],
                                  rmap['ra'],rmap['dec'],1.0/3600.0)
        
        #print('Found {} tractor catalogs in redMaPPer catalog!'.format(len(m1)))
        cat1 = Table(tractor[m1])
        if cat is None:
            cat = cat1
        else:
            cat = vstack((cat,cat1),join_type='exact')
    # Make sure to rm tractor.fits before running because there is no clobber=True
    cat.write(out_dir+'tractor.fits')   
#    fits.writeto(out_dir+'tractor.fits',cat,clobber=True)
    
if __name__ == '__main__':
    main()
