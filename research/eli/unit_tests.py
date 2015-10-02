import numpy as np
import jem_utilities as jem

# Test mag
vec1 = np.array([1.0,1.0,1.0])
m = jem.mag(vec1)
print "Answer should be %.8f" % (np.sqrt(3))
print "jem.mag returns  %.8f" % (m)
