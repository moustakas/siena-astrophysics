#! /usr/bin/env python

from pylab import *
import matplotlib.pyplot as plt
import numpy as np
import sextutils as se
from scipy import optimize


def fitfunc(p,x):
    ret = p[1]*x + p[0]
    return ret


def errfunc(p,x,y,yerr):
    ret = (((fitfunc(p,x)-y)**2)/yerr**2).sum()
    return ret



def bcg_plot(clustername):
    datapath = '/home/obsastro1/siena-astrophysics/summer2013/clash/sbprofiles/'
    savepath = '/home/obsastro1/clash/'+clustername+'_sbprofile.png'

    bcg_file = datapath+clustername+ '_sbprofiles.txt'
    cat = se.se_catalog(bcg_file)

    for i in range(len(cat.sma)):
         x = cat.sma[i]
         y = cat.f160w[i]
         plt.plot(x,y,'bo')
         
    ax = plt.gca()
    ax.set_ylim(ax.get_ylim()[::-1])
    
    plt.title('Light Profile')
    plt.xlabel('SMA')
    plt.ylabel('F160W')
    plt.savefig(savepath)

    plt.show()



def best_fit():
    yerr = np.ones(len(x))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(x,y,yerr= yerr,fmt='o')

    params_starting_vals = [1.0,1.0]
    params_final_vals = optimize.fmin(errfunc,params_starting_vals[:],args=(x,y,yerr),full_output=True)


    
if __name__ == '__main__':
    clustername = 'a1423'
    bcg_plot(clustername)
