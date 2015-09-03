#!/usr/bin/env python

import numpy as np
from astropy.io import ascii
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FixedLocator

def make_sedplot(cat,models,cluster):

    filtername = ['F435W','F606W','F814W','F105W','F125W','F140W','F160W','IRAC_CH1','IRAC_CH2']
    filterwave = [0.43170,0.59177,0.80598,1.0552,1.2486,1.3923,1.5369,3.551,4.496]

    flux = np.array([cat['f435w_flux'],cat['f606w_flux'],cat['f814w_flux'],
                     cat['f105w_flux'],cat['f125w_flux'],cat['f140w_flux'],
                     cat['f160w_flux'],cat['irac_ch1_flux'],cat['irac_ch2_flux']])
    ferr = np.array([cat['f435w_fluxerr'],cat['f606w_fluxerr'],cat['f814w_fluxerr'],
                     cat['f105w_fluxerr'],cat['f125w_fluxerr'],cat['f140w_fluxerr'],
                     cat['f160w_fluxerr'],cat['irac_ch1_fluxerr'],cat['irac_ch2_fluxerr']])
    galaxy = cat['ID']

    fig, ax = plt.subplots()
    plt.rc("font",size=16)
    plt.plot(models['wave']/1E4,models['flux']*10**(0.4*(48.6+31.4)),color='grey') # nJy
    plt.errorbar(filterwave,flux,yerr=ferr,fmt='gs',markersize=15)
    #plt.plot(filterwave,flux,'gs',markersize=15)
    plt.axis([0.2,5.0,-10,np.max(flux)*1.3])
    plt.xlabel(r'Wavelength ($\mu$m)')
    plt.ylabel('Flux (nJy)')
    #plt.loglog()
    #ax.set_xscale('log')
    ax.text(0.05,0.9,cat['ID'],transform=ax.transAxes,fontsize=15)
    #ax.get_xaxis().set_major_locator(FixedLocator(['1.0','3.0','5.0']))
    #ax.get_xaxis().set_major_formatter(ScalarFormatter())
    #ax.set_xticks(['1.0','3.0','5.0'])
    #ax.set_xticklabels(['1.0','3.0','5.0'])

    plt.savefig('seds/'+cluster+'_'+galaxy+'_sed.png')

if __name__ == '__main__':

    topdir = '/Users/ioannis/research/people/ben/'
    for cluster in ['a2744']:
        cat = ascii.read(topdir+cluster+'_hiz.cat',format='sextractor')
        # read the spectral models
        hdu = fits.open(topdir+cluster+'_hiz_models.fits')
        models = hdu[1].data # assuming the first extension is a table
        for ii in range(0,len(cat)):
            make_sedplot(cat[ii],models[ii],cluster)
