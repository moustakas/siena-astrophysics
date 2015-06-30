#! /usr/bin/env python

import os
import numpy as np
from astropy.io import fits
import seaborn as sns
import matplotlib.pyplot as plt


def main():
    """This script simply takes input cluster richness values and outputs halo mass.

    Also outputs a figure of M(BCG) vs M(Halo).

"""
    in_dir = os.getenv('HOME')+'/redmapper/'
    
    isedfit = fits.getdata(in_dir+'decals-isedfit.fits')
    redmapper = fits.getdata(in_dir+'decals-redmapper.fits')
    
    redshift = isedfit['z']
    mbcg = isedfit['mstar_avg']
    lam = redmapper['lambda_chisq']
    mhalo = np.log10((10.0**14.0)*np.exp(1.72+1.08*np.log(lam/60.0)))

    #M(halo) vs M(BCG)
    sns.set(style='white',font_scale=1.5)
    sns.set_context("poster")
    fig = plt.figure(figsize=(8,6))
    sns.kdeplot(mhalo, mbcg)
    #plt.scatter(mhalo, mbcg, c='b')
    plt.xlabel('$log_{10}\ (M(Halo)/M_{sun})$')
    plt.ylim(10,13)
    plt.xlim(13.5,14.8)
    plt.ylabel('$log_{10}\ (M(BCG))$')
    fig.subplots_adjust(bottom=0.2, left=0.15)
    plt.savefig(in_dir+'mbcg_mhalo.jpg',clobber=True)
    plt.close()
    
    #M(BCG) vs redshift
    sns.set_style('white')
    lo = np.where(mhalo<13.9)
    #mid = np.where((mhalo>13.9)*1 and (mhalo<14.5)*1)
    hi = np.where(mhalo>14.5)

    fig2 = plt.figure(figsize=(8,6))
    ax = sns.kdeplot(redshift[lo], mbcg[lo], cmap="Purples", legend=True)
    #middle = sns.kdeplot(redshift[mid],mbcg[mid], cmap="Reds", legend=True, shade_lowest=False, shade=True)
    ax = sns.kdeplot(redshift[hi],mbcg[hi], cmap="Blues", legend=True, ax=None)
    plt.xlabel('redshift')
    plt.ylabel('$log_{10}\ (M(BCG))$')
    plt.ylim(10,13)
    plt.xlim(0,0.8)
    fig.subplots_adjust(bottom=0.2, left=0.15)
    #plt.legend(['low halo', 'mid halo', 'high halo'])
    plt.savefig(in_dir+'3_halo.jpg',clobber=True)

    
if __name__ == '__main__':
    main()
