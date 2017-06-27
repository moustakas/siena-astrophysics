import numpy as np
a=np.random.random(1000000)
b=sum(a)
print b
np.savetxt('elicluster.txt',b)
