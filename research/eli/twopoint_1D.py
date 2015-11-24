import numpy as np
from operator import add
import matplotlib.pylab as plt
import math
from matplotlib.colors import LogNorm
import matplotlib as mpl
DD=np.loadtxt('DD.dat',dtype='float',unpack=True,delimiter=',')
DR=np.loadtxt('DR.dat',dtype='float',unpack=True,delimiter=',')
RR=np.loadtxt('RR.dat',dtype='float',unpack=True,delimiter=',')

ndata=DD[0][0]
nrand=RR[0][0]

x = DD[1][1:]
DD = DD[3][1:]
DR = DR[3][1:]
RR = RR[3][1:]

DD /=(ndata**2-ndata)/2.
DR /=(nrand*ndata)/1.
RR /=(nrand**2-nrand)/2.
w = (DD - 2*DR + RR)/RR

print w

plt.figure(figsize=(8,8))
plt.plot(x,w,'ko')
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title(r'$\xi$')

plt.tight_layout()

plt.show()
