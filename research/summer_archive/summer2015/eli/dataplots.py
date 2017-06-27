import numpy as np
import math
import sys
filename=sys.argv[1]
data=np.loadtxt(filename)
x=data[:,0]
y=data[:,1]
import matplotlib.pyplot as plt
plt.plot(x,y,'go')
plt.xlabel('x Coordinate')
plt.ylabel('y Coordinate')
plt.title('Scatter Plot for Data Set 1')
plt.show()


    
    
