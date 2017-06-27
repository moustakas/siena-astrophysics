import matplotlib.pylab as plt
import numpy as np

import sys

infilename = sys.argv[1]

vals = np.loadtxt(infilename,unpack=True)

ra = vals[0]
dec = vals[1]
redshift = vals[2]

plt.figure(figsize=(10,4))

plt.subplot(1,2,1)
plt.plot(ra,dec,'k.',markersize=0.5)

plt.subplot(1,2,2)
plt.hist(redshift,bins=50)

plt.show()
