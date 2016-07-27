import numpy as np
from operator import add
import matplotlib.pylab as plt
import math
from matplotlib.colors import LogNorm
import matplotlib as mpl
import sys

tag = ""
if len(sys.argv)>1:
    tag = sys.argv[1]

DD=np.loadtxt(tag+'DD.dat',dtype='float',unpack=True,delimiter=',')
DR=np.loadtxt(tag+'DR.dat',dtype='float',unpack=True,delimiter=',')
RR=np.loadtxt(tag+'RR.dat',dtype='float',unpack=True,delimiter=',')

ndata=DD[0][0]
nrand=RR[0][0]

sum_wd = DD[2][0]
sum_wd2 = DD[5][0]
sum_wr = RR[2][0]
sum_wr2 = RR[5][0]

print sum_wd
print sum_wd2 
print sum_wr 
print sum_wr2 
print

norm_dd = 0.5*(sum_wd*sum_wd - sum_wd2) 
norm_rr = 0.5*(sum_wr*sum_wr - sum_wr2)
#norm_dd = sum_wd*sum_wd
#norm_rr = sum_wr*sum_wr
norm_dr = sum_wd*sum_wr

print norm_dd
print norm_dr
print norm_rr

x = DD[1][1:]
DD = DD[3][1:]
DR = DR[3][1:]
RR = RR[3][1:]

print ndata,nrand

#DD /= ((ndata**2-ndata)/2.)
#DR /= ((nrand*ndata)/1.)
#RR /= ((nrand**2-nrand)/2.)
DD /= norm_dd
DR /= norm_dr
RR /= norm_rr
w = (DD - 2*DR + RR)/RR
#w = (DD - 4*DR + RR)/(2*RR)

print DD[-4:]
print DR[-4:]
print RR[-4:]

print w
print ndata
print nrand

plt.figure(figsize=(8,8))
#plt.plot(x,w,'ko',label='Siena')
#plt.ylabel(r'$\xi$',fontsize=24)
plt.plot(x,w*x*x,'ko',label='Siena')
plt.ylabel(r'$\xi r^2$ (Mpc$^2$)',fontsize=24)

plt.xlabel(r'Comoving separation (h$^{-1}$Mpc)',fontsize=18)
plt.title(r'$\xi$')
#plt.xlim(50,180)
#plt.ylim(-0.01,0.10)

# Overlay other results. 
x,y,yerr = np.loadtxt('anderson_results.dat',unpack=True)
#plt.errorbar(x,y,yerr=yerr,fmt='ro',label='Anderson')
plt.errorbar(x,y*x*x,yerr=yerr*x*x,fmt='ro',label='Anderson')

plt.legend()

plt.tight_layout()

plt.show()
