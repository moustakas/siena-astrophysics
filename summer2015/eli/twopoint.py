import numpy as np
from operator import add
import matplotlib.pylab as plt

DD=np.loadtxt('DDtest2d.txt',dtype='float')
DR=np.loadtxt('DRtest2d.txt',dtype='float')
RR=np.loadtxt('RRtest2d.txt',dtype='float')


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

ndata=5000
nrand=10000

#ndata=2
#nrand=2

DD /=(ndata**2-ndata)/2.
DR /=(nrand*ndata)/1.
RR /=(nrand**2-nrand)/2.
theta = (DD - 2*DR + RR)/RR
plt.figure(figsize=(10,10))


#extent=
#plot=plt.imshow(theta)
extent=[-150,150,-150,150]

plt.subplot(2,2,1)
a=plt.imshow(DD,extent=extent)
plt.xlabel('Rperp (Mpc)')
plt.ylabel('Rpara (Mpc)')
plt.title('DD')

plt.subplot(2,2,2)
a=plt.imshow(RR,extent=extent)
plt.xlabel('Rperp (Mpc)')
plt.ylabel('Rpara (Mpc)')
plt.title('RR')

plt.subplot(2,2,3)
a=plt.imshow(DR,extent=extent)
plt.xlabel('Rperp (Mpc)')
plt.ylabel('Rpara (Mpc)')
plt.title('DR')

plt.subplot(2,2,4)
a=plt.imshow(theta,extent=extent)
plt.colorbar(a)
plt.xlabel('Rperp (Mpc)')
plt.ylabel('Rpara (Mpc)')
plt.title('Theta')

#plt.ylim(0,1)
#plt.xlim(0,30)
#plt.xlabel('Distance (Mpc)')
#plt.ylabel('Theta')
#plt.title('DR10 Correlation Estimator with 20000 Galaxies')

plt.tight_layout()

plt.show()


