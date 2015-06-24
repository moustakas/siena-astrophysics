import numpy as np
from operator import add
import matplotlib.pylab as plt
import math
from matplotlib.colors import LogNorm
import matplotlib as mpl
DD=np.loadtxt('DDtest2d4.txt',dtype='float')
DR=np.loadtxt('DRtest2d4.txt',dtype='float')
RR=np.loadtxt('RRtest2d4.txt',dtype='float')


#DDvals=DDvals.transpose()
#DRvals=DRvals.transpose()
#RRvals=RRvals.transpose()

#print DDvals
#exit()

#per = DDvals[0]
#par = DDvals[1]
#DD=DDvals[2]
#DR=DRvals[2]
#RR=RRvals[2]

#print sum(DD)
#print sum(DR)
#print sum(RR)

DD = DD.transpose()
RR = RR.transpose()
DR = DR.transpose()

ndata=150000
nrand=200000


DD /=(ndata**2-ndata)/2.
DR /=(nrand*ndata)/1.
RR /=(nrand**2-nrand)/2.
theta = (DD - 2*DR + RR)/RR

#R^2 WEIGHTING

nbins=200
rangeval=300

# Correct for little h
#rangeval *= 0.7

#R Values
'''
for i in range(nbins):
    for j in range(nbins):
        r2=((nbins/2)-i)**2 + (j-(nbins/2))**2
        theta[i][j] *= r2
'''


plt.figure(figsize=(10,10))


#extent=
#plot=plt.imshow(theta)
extent=[-rangeval,rangeval,-rangeval,rangeval]

plt.subplot(2,2,1)
a=plt.imshow(DD,extent=extent)
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title('DD')

plt.subplot(2,2,2)
b=plt.imshow(RR,extent=extent)
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title('RR')

plt.subplot(2,2,3)
c=plt.imshow(DR,extent=extent)
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title('DR')

### Mirror Over the X-Axis #### 
newtheta= np.zeros((nbins,nbins))
newtheta += theta
for i in range(0,nbins):
    newtheta[i] += theta[(nbins-1)-i]

plt.subplot(2,2,4)
d=plt.imshow(newtheta,extent=extent,norm=mpl.colors.LogNorm(vmin=0.0001,vmax=.2))
plt.colorbar(d)
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title(r'$\xi$')

#plt.ylim(0,1)
#plt.xlim(0,30)
#plt.xlabel('Distance (Mpc)')
#plt.ylabel('Theta')
#plt.title('DR10 Correlation Estimator with 20000 Galaxies')

plt.tight_layout()

plt.show()


