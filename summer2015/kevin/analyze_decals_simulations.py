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
from astropy.io import fits
import seaborn as sns

# Global variables.
scratch_dir = '/global/work/decam/scratch/'
fake_decals_dir = os.getenv('FAKE_DECALS_DIR')

logging.basicConfig(format='%(message)s',level=logging.INFO,stream=sys.stdout)
log = logging.getLogger('decals_simulations')

def flux_graph():
    plt.plot(true_flux,magnitude_difference,'ko',markersize=3)
    plt.xlabel('r (AB mag)')
    plt.ylabel('$\Delta$m (Tractor minus Input, AB mag)')
    plt.ylim([-4,4])

def size_graph():
    size = cat['DISK_R50'][m2]
    plt.plot(size,magnitude_difference,'bo',markersize=3)
    plt.xlabel('$r_{50}$ (arcsec)')
    plt.ylim([-4,4])

def ellipticity_graph():
    ellipticity = cat['DISK_BA'][m2]
    plt.plot(ellipticity,magnitude_difference,'ro',markersize=3)
    plt.xlabel('b/a')
    plt.xlim([0.2,1.0])
    plt.ylim([-4,4])

def ra_position_graph():
    ra_position = cat['ra'][m2]
    plt.plot(ra_position,magnitude_difference,'co',markersize=3)
    plt.xlabel('ra')
    plt.ylim([-4,4])

def dec_position_graph():
    dec_position = cat['dec'][m2]   
    plt.plot(dec_position,magnitude_difference,'yo',markersize=3)
    plt.xlabel('dec')
    plt.ylim([-4,4])

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

    # Make some plots!
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
    tractor_flux = trac['decam_flux'][m1,2]
    true_flux = cat['r'][m2]
    magnitude_difference = -2.5*np.log10(tractor_flux)+22.5-true_flux

    fig = plt.figure(figsize=(8,5))

    plt.subplot(1,3,1)
    flux_graph()

    plt.subplot(1,3,2)
    size_graph()

    plt.subplot(1,3,3)
    ellipticity_graph()

    #plt.subplot(1,5,4)
    #ra_position_graph()
   
    #plt.subplot(1,5,5,sharey=True)
    #dec_position_graph()

    #fig,axes = plt.subplots(ncols=3, sharex=False, sharey=True)

    fig.subplots_adjust(bottom=0.18)
    qafile = os.path.join(fake_decals_dir,'qa_'+brickname+'_resid.png')
    print('Writing {}'.format(qafile))
    plt.savefig(qafile,clobber=True)
    #plt.savefig('/home/desi3/kevin'+'/2428p117_graphs',clobber=True)
    #plt.show()
    
if __name__ == "__main__":
    main()
