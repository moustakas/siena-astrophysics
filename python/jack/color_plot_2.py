#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sextutils as se
import os
import pyfits
import pylab

def makeplot(cat,outfile):
    plt.figure()
    print outfile, len(cat.ra)
#   outpath = '/home/obsastro2/siena-astrophysics/summer2013/clash/'
    for ii in range(len(cat.cluster)):
        clustername = cat.shortname[ii]
#        brightplot = '/home/obsastro2/siena-astrophysics/summer2013/clash/sbprofiles/'+clustername+'_sbprofiles.txt'
        photometry = '/home/obsastro2/clash/apphot/'+clustername+'_bcg_apphot.fits.gz'
        print photometry
        if os.path.isfile(photometry):
            hdu = pyfits.open(photometry)
            data = hdu[1].data
            hdu.close()
            radius = data.field('radius')[0]
            mu = data.field('mu')[0]
            color = mu[:,9]-mu[:,12]
            print clustername, max(color)
            plt.plot(radius,color,'-o', label = clustername)
#           plt.errorbar(cat1.sma,cat1.f110w-cat1.f160w,yerr=cat1.f160w_err,fmt='-o')
    ax = plt.gca()
    ax.set_xscale('log')
    pylab.xlim([.1,10000])
    pylab.ylim([-.1,.8])
    plt.xlabel('Radius (kpc)')
    plt.ylabel('Near-Infrared Color (F110W-F160W)')
    plt.title('Color Profile')
    plt.legend()
    plt.savefig(outfile)
#   plt.show()


if __name__ == "__main__":
    
#   all clusters
    filepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-sample.txt'
    outfile = '/home/obsastro2/siena-astrophysics/summer2013/clash/allclustersplot_test.png'
    cat = se.se_catalog(filepath)
    makeplot(cat,outfile)

    
#   merging galaxies
    mfilepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-merging.txt'
    outfile = '/home/obsastro2/siena-astrophysics/summer2013/clash/mergingclustersplot_test.png'
    mcat = se.se_catalog(mfilepath)
    makeplot(mcat,outfile)


#   star forming filaments
    sfilepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-star_forming.txt'
    outfile = '/home/obsastro2/siena-astrophysics/summer2013/clash/starclustersplot_test.png'
    scat = se.se_catalog(sfilepath)
    makeplot(scat,outfile)

    
#   tidal streams and dust lanes
    cdfilepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-tidal.txt'
    outfile = '/home/obsastro2/siena-astrophysics/summer2013/clash/cdclustersplot_test.png'
    tcat = se.se_catalog(cdfilepath)
    makeplot(tcat,outfile)


#   other galaxies
    ofilepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-other.txt'
    outfile = '/home/obsastro2/siena-astrophysics/summer2013/clash/oclustersplot_test.png'
    ocat = se.se_catalog(ofilepath)
    makeplot(ocat,outfile)
