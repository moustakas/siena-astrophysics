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

cat = ascii.read(codedir+'z9_candidates.cat',format='sextractor')
#print(cat['bpz'], ised['Z'])

mstar = ised['MSTAR_50']#-np.log10(cat['mu'])
mstarerr = ised['MSTAR_ERR']
sfr = ised['SFR_50']#-np.log10(cat['mu'])
sfrerr = ised['SFR_ERR']

fig = plt.figure(figsize=(8,6))
plt.errorbar(mstar,sfr,xerr=mstarerr,yerr=sfrerr,fmt='s',
             ecolor=col[2],capthick=2, markersize=10)
plt.xlabel('$\log_{10}$ (Stellar Mass) $(M_{\odot})$')
plt.ylabel('$\log_{10}$ (Star Formation Rate) $(M_{\odot}\ yr^{-1})$')
#plt.xlim((6.2,9.5))
#plt.ylim((-1.7,1.0))

#plt.text(0.95,0.9,'$<log_{10}$'+' $\sigma>$ = '+
#         '{:.3f}$\pm${:.3f} km/s'.format(gauss.mean.value,gauss.stddev.value),
#         horizontalalignment='right',color='black',
#         transform=plt.gca().transAxes, fontsize=18)
fig.subplots_adjust(left=0.15,bottom=0.15)
fig.savefig(topdir+'sfr_vs_mstar.pdf')
