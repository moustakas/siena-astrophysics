import numpy as np
import math
import sys
filename=sys.argv[1]
data=np.loadtxt(filename)
x=data[:,0]
y=data[:,1]

d = []

for i in range(0,len(data)):     #Cycling through each data point
    stepx=x[i]
    stepy=y[i]
    for n in range(0,len(data)):   #Finding distances between a data point and every other data point
        val=[math.sqrt(((stepx-x[n])**2)+((stepy-y[n])**2))]
        d.append(val[0])
d = np.array(d)

import matplotlib.pyplot as plt
plt.hist(d,50,facecolor='b')
plt.ylabel('Frequency')
plt.xlabel('Distance Values')
plt.title('Data Set 0 Distance Distribution')
plt.show()
