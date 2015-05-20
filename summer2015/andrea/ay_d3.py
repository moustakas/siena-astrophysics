#!/usr/bin/env python


#Number 1 = Gaussian Naive Bayes classification

import numpy as np
from sklearn.naive_bayes import GaussianNB

X=np.random.random((100,2)) #creates 100 points in 2D
y=(X[:,0]+X[:,1]>1).astype(int) #simple division

gnb=GaussianNB()
gnb.fit(X,y)
y_pred=gnb.predict(X)

