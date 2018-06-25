#! /usr/bin/env python

#SVM implementation for both classification and regression tasks.
#Good for discriminative classification 
#At the cost of high contamination level

import numpy as np
from sklearn.svm import LinearSVC

X = np.random.random((100,2))
y = (X[:,0] + X[:,1]>1).astype(int)

model = LinearSVC()
model.fit(X,y)
y_pred = model.predict(X)
