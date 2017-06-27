import numpy as np
data0=np.loadtxt("data1.dat")
column1=data0[:,0]
column2=data0[:,1]
import matplotlib.pyplot as plt
plt.plot(column1,column2,'o')
plt.show()
