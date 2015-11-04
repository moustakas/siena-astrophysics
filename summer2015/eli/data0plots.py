import numpy as np
import math
data0=np.loadtxt("data0.dat")
x=data0[:,0]
y=data0[:,1]
import matplotlib.pyplot as plt
plt.plot(x,y,'o')
plt.show()
print len(x)
#for i in range(0,len(data0)+1):
  #  d=math.sqrt(((x[i]-x[i+1])**2)+((y[i]-y[i+1])**2))
   # print d
    
    
    
