#!/usr/bin/env python

import numpy as np
from glob import glob

from astropy.io import ascii
from astropy.io import fits

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FixedLocator

topdir = '/Users/ioannis/research/people/ben/'

def make_sedplot(cat,sed):

    gal = cat['name']
    lgal = gal.replace('_',' ')
    lgal = lgal.replace(' parr','P')
    print(lgal)
    cluster = cat['cluster']

    filtername = np.array(['F435W','F606W','F814W','F105W','F125W','F140W','F160W','Ks','IRAC_CH1','IRAC_CH2'])
    filterwave = np.array([0.43170,0.59177,0.80598,1.0552,1.2486,1.3923,1.5369,2.14208,3.551,4.496])

    factor = 10**(0.4*31.4)
    flux = sed['maggies']*factor
    ivar = sed['ivarmaggies']/factor**2

    sedwave = sed['wave']/1E4      # micron
    sedflux = sed['flux']*1E23*1E9 # nJy

    fig, ax = plt.subplots()
    plt.rc('font',size=18)
    plt.plot(sedwave,sedflux,color='grey')

    good = np.where((flux>0)*(ivar>0)*1)[0]
    lim = np.where((flux<=0)*(ivar>0)*1)[0]

    ferr = 1.0/np.sqrt(ivar[good])
    fluxlim = 1.0/np.sqrt(ivar[lim])
    #ferr = 1.0/np.sqrt(ivar+((ivar==0)*1))*((ivar!=0)*1)*factor

    mnwave = 0.2
    mxwave = 5.0
    mnflux = -1
    mxflux = np.max(fluxlim)*1.3
    if np.max(flux)>mxflux:
        mxflux = np.max(flux)*1.3
    if gal=='M0416_8958':
        mxflux = 60.0
    if gal=='M0717_11385':
        mxflux = 45.0
    
    plt.errorbar(filterwave[good],flux[good],yerr=ferr,ecolor='dodgerblue',
                 color='dodgerblue',fmt='s',markersize=18,elinewidth=3)
    plt.plot(filterwave[lim],fluxlim,'gv',markersize=18,mfc='none',mew=3)
    #plt.plot(filterwave,flux,'gs',markersize=15)
    plt.axis([mnwave,mxwave,mnflux,mxflux])
    #plt.axis([mnwave,mxwave,0,np.max(sedflux[(sedwave>mnwave)*(sedwave<mxwave)])*1.3])
    #plt.xlabel(r'Observed-Frame Wavelength ($\mu$m)')
    plt.xlabel(r'Observed-Frame Wavelength (micrometers)')
    plt.ylabel('Flux Density (nano-Jansky)')
    #plt.yscale('log')
    #plt.loglog()
    ax.text(0.05,0.88,lgal,transform=ax.transAxes,fontsize=15)
    ax.text(0.05,0.82,'z={:.2f}'.format(cat['bpz']),transform=ax.transAxes,fontsize=15)
    #ax.get_xaxis().set_major_locator(FixedLocator(['1.0','3.0','5.0']))
    #ax.get_xaxis().set_major_formatter(ScalarFormatter())
    #ax.set_xticks(['1.0','3.0','5.0'])
    #ax.set_xticklabels(['1.0','3.0','5.0'])
    fig.tight_layout()

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)

    plt.savefig(topdir+'seds/'+cluster+'_'+gal+'_sed.png')

if __name__ == '__main__':

    # Read the master catalog
    codedir = '/Users/ioannis/repos/git/siena-astrophysics/research/hff/'
    cat = ascii.read(codedir+'z9_candidates.cat',format='sextractor')
    seds = fits.getdata(topdir+'isedfit/ben_fsps_v2.4_miles_chab_none_sfhgrid01_seds.fits.gz', 1)
    
    for ii, obj in enumerate(cat):
        make_sedplot(obj,seds[ii])
