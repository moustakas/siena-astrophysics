#!/usr/bin/pyhton

# Import tools
import numpy as np
import matplotlib.pyplot as plt

# Define path to data file
path='/home/desi1/repos/CUTE/CUTE/test/'

# Load data file
data=np.loadtxt(path+'corr_full_pm.dat')

# Define empty matrices
matrix0=np.zeros(len(data))
matrix1=np.zeros(len(data))
matrix2=np.zeros(len(data))
matrix3=np.zeros(len(data))
matrix4=np.zeros(len(data))
matrix5=np.zeros(len(data))

# Fill matrices with data
for i in range(0,len(data)):
    matrix0[i]=data[i][0]
    matrix1[i]=data[i][1]
    matrix2[i]=data[i][2]
    matrix3[i]=data[i][3]
    matrix4[i]=data[i][4]
    matrix5[i]=data[i][5]

# Plot data on loglog axes and show
plt.loglog(matrix0,matrix1,'bo')
plt.show()
