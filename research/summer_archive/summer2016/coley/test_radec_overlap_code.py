#!/usr/bin/env python

from legacyanalysis.decals_sim import no_overlapping_radec
import matplotlib.pyplot as plt
import numpy as np

bounds = (120.0, 120.25, 30.0, 30.25)
nobj = 100
ra = np.random.uniform(bounds[0], bounds[1], nobj)
dec = np.random.uniform(bounds[2], bounds[3], nobj)

ra, dec = no_overlapping_radec(ra, dec,  bounds, random_state=None, dist=30.0/3600.0) 

#print(ra)
#print(dec)

print(len(ra))

plt.figure()
plt.scatter(ra, dec, s=5, facecolors='none', edgecolors='r')
plt.savefig('testoverlap.png')

