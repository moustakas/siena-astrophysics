#!/usr/bin/env python

"""Label the coadds.
"""

from __future__ import division, print_function

import os
import sys
import logging
import argparse
import numpy as np

from astropy.io import fits
from PIL import Image, ImageDraw

# Global variables.
scratch_dir = '/home/desi3/scratch/'
#scratch_dir = '/global/work/decam/scratch/'
fake_decals_dir = os.getenv('FAKE_DECALS_DIR')

logging.basicConfig(format='%(message)s',level=logging.INFO,stream=sys.stdout)
log = logging.getLogger('decals_simulations')

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='DECaLS simulations.')
    parser.add_argument('-b', '--brick', type=str, default='2428p117', metavar='', 
                        help='process this brick (required input)')

    args = parser.parse_args()
    if args.brick is None:
        parser.print_help()
        sys.exit(1)

    brickname = args.brick
    log.info('Analyzing brick {}'.format(brickname))

    # Read the prior parameters
    priorsfile = os.path.join(fake_decals_dir,'priors_'+brickname+'.fits')
    log.info('Reading {}'.format(priorsfile))
    cat = fits.getdata(priorsfile,1)
    nobj = len(cat)

    #rad = 3*cat['r50_1']/0.262 # half-light radius [pixels]
    rad = np.ones(nobj)*15.0

    imfile = os.path.join(scratch_dir,'coadd',brickname[:3],brickname,'decals-'+brickname+'-image.jpg')
    im = Image.open(imfile)
    sz = im.size
    draw = ImageDraw.Draw(im)
    for ii in range(nobj):
        draw.ellipse((cat['X'][ii]-rad[ii], sz[1]-cat['Y'][ii]-rad[ii],
                      cat['X'][ii]+rad[ii], sz[1]-cat['Y'][ii]+rad[ii]))
    im.save(os.path.join(fake_decals_dir,'qa_'+brickname+'_coadd.png'))
    im.close()

if __name__ == "__main__":
    main()
