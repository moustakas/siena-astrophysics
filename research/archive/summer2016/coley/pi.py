
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
import matplotlib.pyplot as plt
import numpy as np

seed = None
rand = np.random.RandomState(seed)
nobj = 100
ra = rand.uniform(0, 0.25, nobj)
dec = rand.uniform(0, 0.25, nobj)

ngood = 0
while ngood<nobj:

    ra = rand.uniform(0, 0.25, nobj)
    dec = rand.uniform(0, 0.25, nobj)
    coord1 = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    rad = 20*u.arcsec
    indx1, indx2, sep, _ = search_around_sky(coord1, coord1, rad)
    npts = len(indx1)
    doubles = []
    for i in range(0,npts-1):
        if indx1[i+1] == indx1[i]:
            doubles.append(indx1[i])
    newra = np.delete(ra,doubles)
    newdec = np.delete(dec,doubles)
    ng = len(newra)
    ngood = ng + ngood
    
print(ngood)

newraa = newra[:nobj]
newdecc = newdec[:nobj]
plt.scatter(newra, newdec, c='orange')
plt.show()
















'''
plt.scatter(ra, dec, c='orange')
plt.show() 
plt.plot(ra[28], dec[28], 'bx', markersize=10)
plt.plot(ra[74], dec[74], 'rx', markersize=10)


indx3 = set(indx1)
print indx3
indx4 = np.where(indx1 != indx2)[0]
print indx4
np.delete(indx3, indx4)
print indx4
print indx5

indx3 = indx1 != indx2
print indx3
while True in indx3: np.delete(True, indx3)
print indx3
set(indx3)

matrix = np.zeros(2,100)
s = 900
for ii  in range(100):
    r = ii+1
    rx = np.random.random(r)
    ry = np.random.random(r)
    x = rx*s
    y = ry*s
    np.append
    for x in matrix and y in matrix:
         if [x,y] == [x,y]:
        
    #np.append
    
    #disx = np.random.rand(n)
    #disy = np.random.rand(n)
    #plt.scatter(disx,disy)

N = 100
randx = np.random.rand(N)
randy = np.random.rand(N)
colors = np.random.rand(N)
s = 900 #arc sec
#area = np.pi * (15 * np.random.rand(N))**2
ra = randx*s
dec = randy*s

#plt.scatter(ra,dec,s=area,c=colors,alpha=0.5 )
plt.scatter(ra,dec,c=colors,alpha=0.9)
plt.show()
'''
