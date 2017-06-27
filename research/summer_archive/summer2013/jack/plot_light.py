#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sextutils as se

def plot_function(clustername):
    datapath = '/home/obsastro2/siena-astrophysics/summer2013/clash/sbprofiles/'
    outpath = '/home/obsastro2/clash/'+clustername+'_sbprofile.png'
    brightplot = datapath + clustername + '_sbprofiles.txt'
    cat = se.se_catalog(brightplot)
    plt.errorbar(cat.sma,cat.f160w,yerr=cat.f160w_err,fmt='-o')
    #for i in range(len(cat.sma)):
        #x = cat.sma[i]
        #y = cat.f160w[i]
        #plt.plot(x,y,'bo')
        #plt.errorbar(x,y,yerr = 
    plt.xlabel('SMA(kpc)')
    plt.ylabel('SB(mag)')
    plt.title('Light Profile')
#flip y-axis
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    plt.show()
    
    
if __name__ == "__main__":
    clustername = 'a1423'
    plot_function(clustername)
