#! usr/bin/env python

#Linear Discriminant Analysis and Quadratic Discriminant Analysis

import numpy as np
from sklearn.lda import LDA
from sklearn.qda import QDA

X = np.random.random((100,2)) #100 pts in 2D
y = (X[:,0] +X[:,1]>1).astype(int)

lda = LDA()
lda.fit(X,y)
y_pred = lda.predict(X)

qda = QDA()
qda.fit(X,y)
y_pred = qda.predict(X)
