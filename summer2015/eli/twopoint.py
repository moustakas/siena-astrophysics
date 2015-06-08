import numpy as np
from operator import add
import matplotlib.pylab as plt

DDvals=np.loadtxt('DDtest2d.txt',dtype='float')
DRvals=np.loadtxt('DRtest2d.txt',dtype='float')
RRvals=np.loadtxt('RRtest2d.txt',dtype='float')

#DR[:,1]*=2
DDvals=DDvals.transpose()
DRvals=DRvals.transpose()
RRvals=RRvals.transpose()

print DDvals
#exit()

per = DDvals[0]
par = DDvals[1]
DD=DDvals[2]
DR=DRvals[2]
RR=RRvals[2]

print sum(DD)
print sum(DR)
print sum(RR)
ndata=5000
nrand=10000
DD /=(ndata**2-ndata)/2
DR /=(nrand*ndata)
RR /=(nrand**2-nrand)/2


#RRdis[-10:-1]=0.00001
#RRdis[-1]=0.00001
#top1=np.subtract(DDdis,DRdis)
#top2=np.add(top1,RRdis)
#vals=np.array((DD[0]))
#theta=np.divide(top2,RRdis)

theta = (DD - 2*DR + RR)/RR
plt.loglog(per,par,'o')
plt.ylim(0,1)
plt.xlim(0,30)
plt.xlabel('Distance (Mpc)')
plt.ylabel('Theta')
plt.title('DR10 Correlation Estimator with 20000 Galaxies')
plt.show()


