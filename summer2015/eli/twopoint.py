import numpy as np
from operator import add
import matplotlib.pylab as plt

DDvals=np.loadtxt('DDtest.txt',dtype='float')
DRvals=np.loadtxt('DRtest.txt',dtype='float')
RRvals=np.loadtxt('RRtest.txt',dtype='float')

#DR[:,1]*=2
DDvals=DDvals.transpose()
DRvals=DRvals.transpose()
RRvals=RRvals.transpose()

print DDvals
#exit()

x = DDvals[0]
DD=DDvals[1]
DR=DRvals[1]
RR=RRvals[1]

print sum(DD)
print sum(DR)
print sum(RR)
ndata=20000
nrand=20000
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
plt.loglog(x,theta,'o')
plt.ylim(0,1)
plt.xlim(0,30)
plt.xlabel('Distance (Mpc)')
plt.ylabel('Theta')
plt.title('DR10 Correlation Estimator with 20000 Galaxies')
plt.show()


