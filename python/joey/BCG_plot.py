#! /usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sextutils as se

def bcg_plot(clustername):
    datapath = '/home/obsastro1/siena-astrophysics/summer2013/clash/sbprofiles/'
    savepath = '/home/obsastro1/clash/'+clustername+'_sbprofile.png'

    bcg_file = datapath+clustername+ '_sbprofiles.txt'
    cat = se.se_catalog(bcg_file)

    for i in range(len(cat.sma)):
         x = cat.sma[i]
         y = cat.f160w[i]
         plt.plot(x,y,'o')

         plt.xlabel('SMA')
         plt.ylabel('F160W')
         plt.savefig(savepath)

    plt.show()
    
if __name__ == '__main__':
    clustername = 'a1423'
    bcg_plot(clustername)
