#!/usr/bin/env python

"""Analyze the output of decals_simulations.
"""

from __future__ import division, print_function

import os
import sys
import logging
import argparse
import numpy as np
import seaborn as sns

import matplotlib.pyplot as plt

from astropy.io import fits
from astrometry.libkd.spherematch import match_radec

# Global variables.
scratch_dir = '/global/work/decam/scratch/'
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

    # Read the Tractor catalog
    trac = fits.getdata(os.path.join(scratch_dir,'tractor',
                                     brickname[:3],'tractor-'+brickname+'.fits'))

    m1, m2, d12 = match_radec(trac['ra'],trac['dec'],cat['ra'],cat['dec'],1.0/3600.0)

    sns.set(style='white',font_scale=1.5)

    # Plot fraction detected
    fig = plt.figure(figsize=(8,6))
    rminmax = np.array([18.0,24.0]) # get this from the priors table!
    binsz = 0.2
    nbin = long((rminmax[1]-rminmax[0])/binsz)
    yall, bins = np.histogram(cat['r'], bins=nbin, range=rminmax)
    ymatch, _ = np.histogram(cat['r'][m2], bins=nbin, range=rminmax)
    cbins = (bins[:-1] + bins[1:]) / 2.0
    plt.plot(cbins,1.0*ymatch/yall,lw=3)
    plt.axhline(y=1.0,lw=2,ls='--',color='gray')
    plt.xlabel('Input r magnitude (AB mag)')
    plt.ylabel('Fraction Detected')
    plt.ylim([0.0,1.1])
    fig.subplots_adjust(bottom=0.15)
    plt.savefig(os.path.join(fake_decals_dir,'qa_'+brickname+'_frac.png'))

    # Residual plots
    


if __name__ == "__main__":
    main()
