import numpy as np
import matplotlib
import matplotlib.pylab as plt

# Generate some fake data
data0 = np.random.random((100,100))
data1 = np.random.random((100,100))

for i in range(0,100):
    for j in range(0,100):
        r = (i-50)**2 + (j-50)**2
        r = np.sqrt(r)
        print r
        if r<40 or r>60:
            data1[i][j]=0



# Write data to files
np.savetxt('data0.dat',data0)
np.savetxt('data1.dat',data1)

# Read it back in
indata0 = np.loadtxt('data0.dat')
indata1 = np.loadtxt('data1.dat')

# Plot
plt.figure(figsize=(16,4))
# Assume that the x-axis and y-axis range from -100,100
extent = [-100,100,-100,100]

ax0 = plt.subplot(1,3,1)
cs0 = ax0.imshow(indata0,extent=extent,interpolation='nearest',origin='lower',cmap=plt.cm.coolwarm,axes=ax0,aspect='auto') 
plt.colorbar(cs0)

ax1 = plt.subplot(1,3,2)
cs1 = ax1.imshow(indata1,extent=extent,interpolation='nearest',origin='lower',cmap=plt.cm.coolwarm,axes=ax1,aspect='auto') 
plt.colorbar(cs1)

ratio = indata1/indata0

ax2 = plt.subplot(1,3,3)
norm = matplotlib.colors.Normalize(vmin=0.0,vmax=20.0)
cs2 = ax2.imshow(ratio,extent=extent,interpolation='nearest',origin='lower',cmap=plt.cm.coolwarm,axes=ax2,aspect='auto',norm=norm)
plt.colorbar(cs2)


plt.show()
