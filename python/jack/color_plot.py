#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sextutils as se
import os

def makeplot(cat,outfile):
    plt.figure()
    print outfile, len(cat.ra)
#   outpath = '/home/obsastro2/siena-astrophysics/summer2013/clash/'
    for ii in range(len(cat.cluster)):
        clustername = cat.shortname[ii]
        brightplot = '/home/obsastro2/siena-astrophysics/summer2013/clash/sbprofiles/'+clustername+'_sbprofiles.txt'
        if os.path.isfile(brightplot):
            cat1 = se.se_catalog(brightplot)
            plt.plot(cat1.sma,cat1.f110w-cat1.f160w,'-o', label = clustername)
#           plt.errorbar(cat1.sma,cat1.f110w-cat1.f160w,yerr=cat1.f160w_err,fmt='-o')
    plt.xlabel('SMA (kpc)')
    plt.ylabel('F110W-F160W')
    plt.title('Color Profile')
    plt.legend()
    plt.savefig(outfile)
#   plt.show()

if __name__ == "__main__":
    
#all clusters
    filepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-sample.txt'
    outfile = '/home/obsastro2/siena-astrophysics/summer2013/clash/allclustersplot.png'
    cat = se.se_catalog(filepath)
    makeplot(cat,outfile)

    
#merging galaxies
    mfilepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-merging.txt'
    outfile = '/home/obsastro2/siena-astrophysics/summer2013/clash/mergingclustersplot.png'
    mcat = se.se_catalog(mfilepath)
    makeplot(mcat,outfile)


#star forming filaments
    sfilepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-star_forming.txt'
    outfile = '/home/obsastro2/siena-astrophysics/summer2013/clash/starclustersplot.png'
    scat = se.se_catalog(sfilepath)
    makeplot(scat,outfile)

    
#tidal streams and dust lanes
    cdfilepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-tidal.txt'
    outfile = '/home/obsastro2/siena-astrophysics/summer2013/clash/cdclustersplot.png'
    tcat = se.se_catalog(cdfilepath)
    makeplot(tcat,outfile)


#other galaxies
    ofilepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-other.txt'
    outfile = '/home/obsastro2/siena-astrophysics/summer2013/clash/oclustersplot.png'
    ocat = se.se_catalog(ofilepath)
    makeplot(ocat,outfile)
