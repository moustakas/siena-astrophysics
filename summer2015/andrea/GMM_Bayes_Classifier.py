#! /usr/bin/env python

#the implimentation of Gaussian Mixture Model as a classifier

import numpy as np
from astroML.classification import GMMBayes

X = np.random.random((100,2))  #100 points in 2D
y = (X[:,0] + X[:,1]>1).astype(int)  #simple division

gmmb = GMMBayes(3)  #three clusters per class
gmmb.fit(X,y)
y_pred = gmmb.predict(X)
