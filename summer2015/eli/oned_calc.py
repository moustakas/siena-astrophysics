import astropy.io 
from astropy.io import fits
import numpy as np
import sys
import random
import math
import scipy
import scipy.spatial
import matplotlib.pylab as plt
from astropy.cosmology import FlatLambdaCDM
import time
import argparse


# Comoving Distances In Mpc
cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
comdist=cosmo.comoving_distance(redshift).value

# Convert spherical to Cartesian Coordinates
x=comdist*np.sin(dec)*np.cos(ra)
y=comdist*np.sin(dec)*np.sin(ra)
z=comdist*np.cos(dec)

coords=np.column_stack((x,y,z))
