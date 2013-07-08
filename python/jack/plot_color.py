#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sextutils as se
import os

if __name__ == "__main__":
    
#all clusters
    filepath = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-sample.txt'
    cat = se.se_catalog(filepath)


    outpath = '/home/obsastro2/siena-astrophysics/summer2013/clash/'
    for ii in range(len(cat.cluster)):
        clustername = cat.shortname[ii]
        brightplot = '/home/obsastro2/siena-astrophysics/summer2013/clash/sbprofiles/'+clustername+'_sbprofiles.txt'
        if os.path.isfile(brightplot):
            cat1 = se.se_catalog(brightplot)
            plt.plot(cat1.sma,cat1.f110w-cat1.f160w, '-o')
            #plt.errorbar(cat1.sma,cat1.f110w-cat1.f160w,yerr=cat1.f160w_err,fmt='-o')
            
            plt.xlabel('SMA (kpc)')
            plt.ylabel('F110W-F160W')
            plt.title('Light Profile')
    plt.show()





