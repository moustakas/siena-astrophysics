import numpy as np

from time import time

# Example of how to time a process
'''
start = time()

for i in xrange(100000):
    1

print "Time to run: %f ms" % (time()-start)
'''


n = 100000000

x = np.random.random(n)

print x[0:10]

xnew = []

# Grab the numbers less than 0.5 the pure python way.
start = time()
for i in range(0,n):
    if x[i]<0.5:
        xnew.append(x[i])
print "Time to run: %f ms" % (time()-start)

print len(xnew)


# Grab the numbers less than 0.5 the numpy/slicing way
start = time()
index0 = x<0.5
xnew = x[index0]
print "Time to run: %f ms" % (time()-start)

print len(xnew)
