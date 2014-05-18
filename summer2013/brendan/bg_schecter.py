from pylab import *
import atpy

lum = arange(0.01,5,0.1)
alpha = [-.5,-.75,-1]
figure()
f=0
while f<3:
    ax=gca()
    ax.set_xscale('log')
    ax=gca()
    ax.set_yscale('log')
    schectlum = (8.*10**-3)*((lum)**(alpha[f]))*(exp(-lum))
    f=f+1
    plot(lum,schectlum)
