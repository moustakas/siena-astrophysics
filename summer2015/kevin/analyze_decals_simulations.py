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

logging.basicConfig(format='%(message)s',level=logging.INFO,stream=sys.stdout)
log = logging.getLogger('decals_simulations')

def get_simdir(brickname=None,objtype=None):
    """Get the simulation directory."""

    # Check for the environment variables we need.
    if 'DECALS_SIM_DIR' in os.environ:
        decals_sim_dir = os.getenv('DECALS_SIM_DIR')
    else:
        log.error('Missing $DECALS_SIM_DIR environment variable')
        sys.exit(1)

    log.info('DECALS_SIM_DIR {}'.format(decals_sim_dir))
    simdir = os.path.join(decals_sim_dir,objtype.lower(),brickname)
    
    return simdir

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='DECaLS simulations.')
    parser.add_argument('-b', '--brick', type=str, default='2428p117', metavar='', 
                        help='process this brick (required input)')
    parser.add_argument('-o', '--objtype', type=str, default='ELG', metavar='', 
                        help='object type (STAR, ELG, LRG, BGS)') 

    args = parser.parse_args()
    if args.brick is None:
        parser.print_help()
        sys.exit(1)

    brickname = args.brick
    objtype = args.objtype.upper()
    lobjtype = objtype.lower()
    log.info('Analyzing objtype {} on brick {}'.format(objtype,brickname))

    decals_sim_dir = get_simdir(brickname,objtype)
    
    # Read the simulated object catalog
    simcatfile = os.path.join(decals_sim_dir,'simcat-'+brickname+
                              '-'+lobjtype+'.fits')
    log.info('Reading {}'.format(simcatfile))
    simcat = fits.getdata(simcatfile,1)

    rminmax = np.array([18.0,25.0]) # get this from the priors table!

    # Read the Tractor catalog
    tractorfile = os.path.join(decals_sim_dir,'tractor-'+
                               brickname+'-'+lobjtype+'.fits')
    log.info('Reading {}'.format(tractorfile))
    tractor = fits.getdata(tractorfile)

    m1, m2, d12 = match_radec(tractor['ra'],tractor['dec'],
                              simcat['ra'],simcat['dec'],1.0/3600.0)

    # Make some plots!
    sns.set(style='white',font_scale=1.6)

    # Residual plots - deal with negative fluxes!!
    rmag_tractor = tractor['decam_flux'][m1,2]
    rmag = simcat['r'][m2]
    deltam = -2.5*np.log10(rmag_tractor)+22.5-rmag

    fig = plt.figure(figsize=(8,5))
    ax = fig.gca()
    plt.ylim([-3,3])
    plt.ylabel('$\Delta$m (Tractor minus Input, AB mag)')
    
    if objtype=='STAR':
        plt.subplot(1,1,1)
        plt.plot(rmag,deltam,'s',markersize=4)
        plt.axhline(y=0.0,lw=2,ls='solid',color='gray')
        plt.xlim(rminmax)
        plt.xlabel('Input r magnitude (AB mag)')
        ax.text(0.05,0.05,objtype,horizontalalignment='left',
                verticalalignment='bottom',transform=ax.transAxes,
                fontsize=16)

    if objtype=='ELG':
        fig = plt.figure(figsize=(8,4))
        plt.subplot(1,3,1)
        plt.plot(rmag,deltam,'s',markersize=3)
        plt.axhline(y=0.0,lw=2,ls='solid',color='gray')
        plt.xlim(rminmax)
        plt.xlabel('r (AB mag)')

        plt.subplot(1,3,2)
        plt.plot(simcat['R50_1'][m2],deltam,'s',markersize=3)
        plt.axhline(y=0.0,lw=2,ls='solid',color='gray')
        plt.xlabel('$r_{50}$ (arcsec)')

        plt.subplot(1,3,3)
        plt.plot(simcat['BA_1'][m2],deltam,'s',markersize=3)
        plt.axhline(y=0.0,lw=2,ls='solid',color='gray')
        plt.xlabel('b/a')
        plt.xlim([0.2,1.0])

    fig.subplots_adjust(bottom=0.18)
    qafile = os.path.join(decals_sim_dir,'qa-'+brickname+'-'+
                          lobjtype+'-deltam.png')
    log.info('Writing {}'.format(qafile))
    plt.savefig(qafile)
    
    # Plot fraction detected
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca()
    binsz = 0.2
    nbin = long((rminmax[1]-rminmax[0])/binsz)
    yall, bins = np.histogram(simcat['r'], bins=nbin, range=rminmax)
    ymatch, _bins = np.histogram(simcat['r'][m2], bins=nbin, range=rminmax)

    good = np.where((np.abs(deltam)<0.3)*1)
    ymatchgood, _binsgood = np.histogram(simcat['r'][m2[good]], bins=nbin, range=rminmax)
    
    cbins = (bins[:-1] + bins[1:]) / 2.0
    plt.step(cbins,1.0*ymatch/yall,lw=3,alpha=0.5,label='All '+objtype+'s')
    plt.step(cbins,1.0*ymatchgood/yall,lw=3,ls='dashed',label='|$\Delta$m|<0.3')
    plt.axhline(y=1.0,lw=2,ls='dotted',color='gray')
    plt.xlabel('Input r magnitude (AB mag)')
    plt.ylabel('Fraction Detected')
    plt.ylim([0.0,1.1])
    plt.legend(loc='lower left')
    #plt.text(0.1,0.1,horizontalalignment='left',verticalalignment='bottom',
    #         transform=ax.transAxes)
    
    fig.subplots_adjust(bottom=0.15)
    qafile = os.path.join(decals_sim_dir,'qa-'+brickname+'-'+lobjtype+'-frac.png')
    log.info('Writing {}'.format(qafile))
    plt.savefig(qafile)

if __name__ == "__main__":
    main()
