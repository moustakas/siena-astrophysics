#! /usr/bin/env python

#A fast K-Neighbors classifier built on a ball tree for fast neighbor searches

import numpy as np
from sklearn.neighbors import KNeighborsClassifier

X = np.random.random((100,2))  #100 pts in 2D
y = (X[:,0] + X[:,1]>1).astype(int)

knc = KNeighborsClassifier(5)  #uses the 5 nearest neighbors
knc.fit(X,y)
y_pred = knc.predict(X)

