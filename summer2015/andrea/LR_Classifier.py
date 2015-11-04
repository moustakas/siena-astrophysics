#! /usr/bin/env python

#Implemetation of logistic regression

import numpy as np
from sklearn.linear_model import LogisticRegression

X = np.random.random((100,2))
y = (X[:,0] + X[:,1]>1).astype(int)

logr = LogisticRegression(penalty = '12')
logr.fit(X,y)
y_pred = logr.predict(X)
