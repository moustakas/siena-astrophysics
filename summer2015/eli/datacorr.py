import numpy as np
import math
import scipy.spatial
'''
r=np.loadtxt('data (1).dat')
ra=r[:,0]
dec=r[:,1]
coords=np.column_stack((ra,dec))
d1=scipy.spatial.distance.pdist(coords)

# Histogram
import matplotlib.pylab as plt
hist=plt.hist(d1,bins=30,range=(0,6))
plt.xlabel('Galactic Distances (Mpc)')
plt.ylabel('Frequency')
plt.title('Galactic Distance Distribution of 5000 Random Mock Galaxies')
plt.show()                
frequ=hist[0]
dist=hist[1]
m=dist[0:-1]
diff=np.diff(dist)/2
mid=m+diff
vals0=np.column_stack((mid,frequ))
np.savetxt('DDtest.txt',vals0)
###########################################

r1=np.loadtxt('random.dat')
ra1=r1[:,0]
dec1=r1[:,1]
coords1=np.column_stack((ra1,dec1))
d2=scipy.spatial.distance.pdist(coords1)

# Histogram
import matplotlib.pylab as plt
hist1=plt.hist(d2,bins=30,range=(0,6))
plt.xlabel('Galactic Distances (Mpc)')
plt.ylabel('Frequency')
plt.title('Galactic Distance Distribution of 5000 Random Mock Galaxies')
plt.show()                
frequ1=hist1[0]
dist1=hist1[1]
m1=dist1[0:-1]
diff1=np.diff(dist1)/2
mid1=m1+diff1
vals1=np.column_stack((mid1,frequ1))
np.savetxt('RRtest.txt',vals1)

###################################################
val=scipy.spatial.distance.cdist(coords,coords1)
vflat = val.flatten()
d3 = vflat

hist2=plt.hist(d3,bins=30,range=(0,6))
plt.xlabel('Galactic Distances (Mpc)')
plt.ylabel('Frequency')
plt.title('DR with 5000 Galaxies')
plt.show()
frequ2=hist2[0]
dist2=hist2[1]
m2=dist2[0:-1]
diff2=np.diff(dist2)/2
mid2=m2+diff2
vals2=np.column_stack((mid2,frequ2))
np.savetxt('DRtest.txt',vals2)

'''
###################################################
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

xperp = DDvals[0]
xpar= DDvals[1]
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
plt.logloger(x,theta,'o')
#plt.ylim(0,1)
#plt.xlim(0,30)
plt.xlabel('Distance (Mpc)')
plt.ylabel('Theta')
plt.title('Correlation Estimator')
plt.show()


