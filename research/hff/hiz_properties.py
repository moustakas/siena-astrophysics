#!/usr/bin/env python

"""
Plot up the physical properties of the sample.
"""

import os
import numpy as np
from glob import glob

from astropy.io import ascii
from astropy.io import fits

import matplotlib
import matplotlib.pyplot as plt

import seaborn as sns
sns.set(style="darkgrid", font_scale=2)
col = sns.color_palette('dark')
matplotlib.rcParams['text.usetex'] = True

topdir = '/Users/ioannis/research/people/ben/isedfit/'
codedir = '/Users/ioannis/repos/git/siena-astrophysics/research/hff/'

isedfile = os.path.join(topdir,'ben_fsps_v2.4_miles_chab_none_sfhgrid01.fits.gz')
ised = fits.getdata(isedfile, 1)

# Read the z=9 candidates
cat = ascii.read(codedir+'z9_candidates.cat',format='sextractor')
#print(cat['bpz'], ised['Z'])

mstar = ised['MSTAR_50']#-np.log10(cat['mu'])
mstarerr = ised['MSTAR_ERR']
sfr = ised['SFR_50']#-np.log10(cat['mu'])
sfrerr = ised['SFR_ERR']

# Read the z=7 candidates
z7cat = ascii.read(codedir+'z7_sfrmass.cat',format='sextractor')

z7mstar = z7cat['MSTAR']
z7sfr = z7cat['SFR']
z7mstarerr = (z7cat['MSTAR_UPERR']+z7cat['MSTAR_LOERR'])/2.0
z7sfrerr = (z7cat['SFR_UPERR']+z7cat['SFR_LOERR'])/2.0

# Figure w/ just z=9

fig, ax = plt.subplots(1, 1, figsize=(8,6))
ax.errorbar(mstar,sfr,xerr=mstarerr,yerr=sfrerr,fmt='s',
             ecolor=col[2],capthick=2, markersize=10, label=r'$z\approx9$')
ax.set_xlabel('$\log_{10}$ (Stellar Mass) $(M_{\odot})$')
ax.set_ylabel('$\log_{10}$ (Star Formation Rate) $(M_{\odot}\ yr^{-1})$')
ax.set_xlim((7.3,9.5))
leg = ax.legend(loc='lower right', frameon=True)
leg.get_frame().set_edgecolor('black')
fig.subplots_adjust(left=0.15,bottom=0.15)
fig.savefig(topdir+'sfr_vs_mstar_z9.pdf')

# Figure w/ z=7

fig, ax = plt.subplots(1, 1, figsize=(8,6))
ax.errorbar(mstar,sfr,xerr=mstarerr,yerr=sfrerr,fmt='s',
             ecolor=col[2],capthick=2, markersize=10, label=r'$z\approx9$')
ax.errorbar(z7mstar,z7sfr,xerr=z7mstarerr,yerr=z7sfrerr,fmt='o',
             ecolor=col[4],capthick=2, markersize=10, label=r'$z\approx7$')
ax.set_xlabel('$\log_{10}$ (Stellar Mass) $(M_{\odot})$')
ax.set_ylabel('$\log_{10}$ (Star Formation Rate) $(M_{\odot}\ yr^{-1})$')
ax.set_xlim((7.3,9.5))
leg = ax.legend(loc='lower right', frameon=True)
leg.get_frame().set_edgecolor('black')
fig.subplots_adjust(left=0.15,bottom=0.15)
fig.savefig(topdir+'sfr_vs_mstar_z7.pdf')
