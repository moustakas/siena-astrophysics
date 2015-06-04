import numpy as np
x=np.array([3,4,65354,64,32,625,57,34])
xin=[]
for i in range (0,len(x)):
    x1=np.binary_repr(x[i])
    xin.append(x1)
xin=np.array(xin)

    
