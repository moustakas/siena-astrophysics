#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from matplotlib.colors import LogNorm

def make_colorcolor(cluster,topdir):
    
    print 'Reading '+topdir+cluster+'_multicolor_nir.cat'
    cat = ascii.read(topdir+cluster+'_multicolor_nir.cat',
                     format='sextractor')

    # find our high-z candidates
    if cluster == 'a2744':
        hiz = ascii.read(topdir+cluster+'_hiz_z9.cat',
                         format='sextractor')
        hiz_f140w_f160w = -2.5*np.log10(hiz['f140w_flux']/hiz['f160w_flux'])
        hiz_f105w_f140w = -2.5*np.log10(hiz['f105w_flux']/hiz['f140w_flux'])

    # now make the plot
    ra = cat['alpha_j2000']
    dec = cat['delta_j2000']
    f105w = cat['hst_wfc3_ir_f105w_mag_bpz']
    f140w = cat['hst_wfc3_ir_f140w_mag_bpz']
    f160w = cat['hst_wfc3_ir_f160w_mag_bpz']

    f105w_err = cat['hst_wfc3_ir_f105w_magerr_bpz']
    f140w_err = cat['hst_wfc3_ir_f140w_magerr_bpz']
    f160w_err = cat['hst_wfc3_ir_f160w_magerr_bpz']
    
    xmin = -3.0
    xmax = 4.0
    ymin = -3.0
    ymax = 5.0
    coeff = [1.0,0.8]
    fitmin = (0.8-coeff[1])/coeff[0]
    xcolor = np.linspace(fitmin,0.6,20)
    ycolor = np.polyval(coeff,xcolor)
    
    plt.rc("font",size=16)
    plt.errorbar(f140w-f160w,f105w-f140w,f140w_err,f140w_err,'bo')
    plt.plot(hiz_f140w_f160w,hiz_f105w_f140w,'gs',markersize=15)
#   plt.plot(f140w-f160w,f105w-f140w,'bo')
    plt.plot(xcolor,ycolor,'r-',[0.6,0.6],
             [np.polyval(coeff,0.6),ymax],'r-',
             [xmin,fitmin],[0.8,0.8],'r-',linewidth=2)
    #plt.hist2d(f140w-f160w,f105w-f140w,bins=50,norm=LogNorm(),
    #range=[[xmin,xmax],[ymin,ymax]])
    #plt.colorbar()
    
    #plt.plot(-2.5*np.log10(21.8/35.4),-2.5*np.log10(1.2/21.8),'gs') # JD1A
    #plt.plot(-2.5*np.log10(24.3/43.8),-2.5*np.log10(2.0/24.3),'gs') # JD1B
    
    plt.axis([xmin,xmax,ymin,ymax])
    plt.xlabel('F140W - F160W')
    plt.ylabel('F105W - F140W')
    plt.savefig(topdir+cluster+'_z9_colorcolor.png')


if __name__ == '__main__':

    topdir = '/Users/ioannis/research/people/ben/'
    for cat in ['a2744']:
    #for cat in ['a2744','m0416']:
        make_colorcolor(cat,topdir)
