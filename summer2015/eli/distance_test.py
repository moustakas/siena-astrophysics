############ DISTANCE TEST ##################
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.io
import matplotlib.pylab as plt
import math

# 10 Galaxies (ra,dec,z)    
galdat=np.array([[221.42511802,6.48006286,0.5891996],[143.26441531,16.56591479,0.68778771],[250.39513305,32.56592786,0.51523036],[192.75572856,-3.06978128,0.56942481],[136.94757689,0.61384908,0.53580755],[126.64582688,13.7770392,0.54264402],[127.16685388,42.06680514,0.53326231],[130.6626883,4.63926095,0.62787837],[242.14836333,12.94153412,0.46952653],[236.24144054,23.61780059,0.53417909]])
ra=galdat[:,0] # Right Ascension Values
dec=galdat[:,1] # Declination Values
z=galdat[:,2] # Redshift

# Comoving Distances
cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
comdist=cosmo.comoving_distance(z)






