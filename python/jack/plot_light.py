#! /usr/bin/env python

#import matplotlib.pyplot as plt
#import numpy as np
#import sextutils as se

def plot_function(x,y):
    rootpath = '/home/obsastro2/siena-astrophysics/'
    inpath = 'summer2013/clash/sbprofiles/'
    outpath = rootpath + 'summer2013/clash/'
    #myplot = se.se_catalog(inpath)
    brightplot = inpath + 'a1423_sbprofiles.txt'
    plot(brightplot.SMA,brightplot.F160w)
    xlabel('SMA(kpc)')
    ylabel('SB(mag)')
    title('Light Profile')

    show()


#if __name__ == "__main__":
            
    #outfile = outpath + str(plot.clustername)
    
