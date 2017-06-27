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
    for ii in range(len(cat.cluster)):
        clustername = cat.shortname[ii]
        photometry = '/home/obsastro2/clash/isedfit/'+clustername+'_fsps_salp_none_sfhgrid01.fits.gz'
        print photometry
        if os.path.isfile(photometry):
            hdu = pyfits.open(photometry)
            data = hdu[1].data
            hdu.close()
            age = data.field('sfrage_50')
            radius = np.log10(data.field('radius'))
            error = data.field('sfrage_err')
            index = np.where(radius<=np.log10(cat.rmax[ii]))
            plt.plot(radius[index],age[index],'-o', label = clustername)
#            plt.errorbar(radius,age,yerr=error,fmt='-o', label = clustername)
#    ax = plt.gca()
#    ax.set_xscale('log')
    pylab.xlim([-.6,2.5])
    pylab.ylim([0,10])
    plt.ylabel('Star Formation Rate Weighted Age (Gyr)')
    plt.xlabel('log10 Radius/semi-major axis (kpc)')
#    plt.title('Star Formation vs Distance')
    plt.legend(loc='lower right')
    fig = plt.gcf()
#    fig.set_size_inches(8,5)
    plt.savefig(outfile)
#   plt.show()


if __name__ == "__main__":
    
#   all clusters
    filepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-sample.txt'
    outfile = '/home/obsastro2/clash/isedfit/allclusters.pdf'
    cat = se.se_catalog(filepath)
    makeplot(cat,outfile)

    
#   merging galaxies
    mfilepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-merging.txt'
    outfile = '/home/obsastro2/clash/isedfit/mergingclusters.pdf'
    mcat = se.se_catalog(mfilepath)
    makeplot(mcat,outfile)


#   star forming filaments
    sfilepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-star_forming.txt'
    outfile = '/home/obsastro2/clash/isedfit/starclusters.pdf'
    scat = se.se_catalog(sfilepath)
    makeplot(scat,outfile)

    
#   tidal streams and dust lanes
    cdfilepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-tidal.txt'
    outfile = '/home/obsastro2/clash/isedfit/cdclusters.pdf'
    tcat = se.se_catalog(cdfilepath)
    makeplot(tcat,outfile)


#   other galaxies
    ofilepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-other.txt'
    outfile = '/home/obsastro2/clash/isedfit/otherclusters.pdf'
    ocat = se.se_catalog(ofilepath)
    makeplot(ocat,outfile)
