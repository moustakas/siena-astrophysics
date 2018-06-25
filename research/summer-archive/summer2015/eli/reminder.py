import numpy as np
testdata=np.loadtxt("test_file.dat.txt")
column1=testdata[:,0]
column2=testdata[:,1]
column3=testdata[:,2]
import matplotlib.pyplot as plt
plt.plot(column1,column2,'o',markersize=0.5)
plt.xlabel('RA (Degrees)')
plt.ylabel('Dec (Degrees)')
plt.title('Reminder Test Data')
plt.show()
