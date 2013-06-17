#!/usr/bin/env python
from pylab import *
import atpy

def schecter():
    f=0
    figure()
    alpha = [-0.5, -0.75, -1.]
    lum =arange(.01, 5., .1)
    ax=gca()
    ax.set_yscale('log')
    ax=gca()
    ax.set_xscale('log')
    while f<3:
        schect = (8.*10**-3)*(lum**alpha[f])*exp(-lum)
        plot(lum, schect)
        f=f+1
