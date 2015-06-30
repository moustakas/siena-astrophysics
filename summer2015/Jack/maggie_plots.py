#! /usr/bin/env python

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os

def main():
    """Document me.

    """
    out_dir = os.getenv('HOME')+'/redmapper/'
    tractor = fits.getdata(out_dir+'tractor.fits',1)

    gg = tractor['decam_flux'][:,1]
    rr = tractor['decam_flux'][:,2]
    zz = tractor['decam_flux'][:,4]

    sdss_gg = tractor['sdss_modelflux'][:,1]
    sdss_rr = tractor['sdss_modelflux'][:,2]
    sdss_zz = tractor['sdss_modelflux'][:,4]

    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
    y_gg = -2.5*np.log10(gg/sdss_gg)
    x_gg = -2.5*np.log10(sdss_gg)+22.5
    y_rr = -2.5*np.log10(rr/sdss_rr)
    x_rr = -2.5*np.log10(sdss_rr)+22.5
    y_zz = -2.5*np.log10(zz/sdss_zz)
    x_zz = -2.5*np.log10(sdss_zz)+22.5

    ax1.scatter(x_gg, y_gg, c="g")
    ax1.set_xlim(13,25)
    ax1.set_ylim(-3,3)
    ax1.axhline(y=0.000,xmin=0,xmax=3,c="black",linewidth=0.5,zorder=0)
    ax2.scatter(x_rr, y_rr, c="r")
    ax2.axhline(y=0.000,xmin=0,xmax=3,c="black",linewidth=0.5,zorder=0)
    ax2.set_ylabel('$\Delta$m (AB mag)')
    ax3.scatter(x_zz, 2*y_zz**2-1,c="b")
    ax3.axhline(y=0.000,xmin=0,xmax=3,c="black",linewidth=0.5,zorder=0)
    ax3.set_xlabel('m_sdss (AB mag)')
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.savefig(out_dir+'maggie_phot.png')
    plt.show()

    
if __name__ == '__main__':
    main()
