#!/usr/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from sklearn import svm, datasets

#iris = datasets.load_iris()
#X = iris.data[:, :2]  # we only take the first two features. We could
#Y = iris.target

filename = '/home/desi2/stars_qsos_sdss.fits'
data = fits.getdata(filename,1)
ug = data['u'] - data['g']
gr = data['g'] - data['r']
objtype = data['class']

colors = np.vstack((ug,gr)).T
labels = (objtype == 'STAR')*1

# we create an instance of SVM and fit out data.
# here we could also put 
#clf = svm.SVC(kernel='poly',degree=4)
clf = svm.SVC(kernel = 'rbf')
#clf = svm.SVC(kernel=my_kernel)
dofit = clf.fit(colors,labels)

# Plot the decision boundary. For that, we will assign a color to each
# point in the mesh [x_min, m_max]x[y_min, y_max].
ugmin = -0.5
ugmax = 3.0
grmin = -0.5
grmax = 1.5

hh = 0.02  # step size in the mesh
xx, yy = np.meshgrid(np.arange(grmin,grmax,hh),np.arange(ugmin,ugmax,hh))
zz = dofit.predict(np.c_[xx.ravel(),yy.ravel()])

# Put the result into a color plot
zz = zz.reshape(xx.shape)

fig = plt.figure()
plt.pcolormesh(xx,yy,zz,cmap=plt.cm.Paired)
plt.scatter(colors[:,1],colors[:,0],c=labels)# ,cmap=plt.cm.Paired)
plt.title('Title')
plt.axis('tight')
plt.xlim([grmin,grmax])
plt.ylim([ugmin,ugmax])
plt.savefig('/home/desi2/sdss_SVM.png')
