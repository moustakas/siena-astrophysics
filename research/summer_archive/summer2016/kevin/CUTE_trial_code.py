#!/usr/bin/pyhton

# Import tools
import numpy as np
import matplotlib.pyplot as plt

# Define path to data file
path='/home/desi1/repos/CUTE/CUTE/test/'

# Load data file
data=np.loadtxt(path+'corr_full_pm.dat')

for ii in range(0,len(data[0])):
    locals()['matrix{0}'.format(ii)]=data[:,ii]

# Plot data on loglog axes and show
plt.loglog(matrix0,matrix1,'bo')
plt.show()


