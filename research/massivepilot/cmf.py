#!/usr/bin/python

import os

import numpy as np
import matplotlib.pyplot as plt

from astropy.io.ascii import read
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import fitsio
import seaborn

# Reading in the RedMaPPer iSEDfit catalog.
rmpath = os.path.join(os.sep, 'global', 'work', 'projects', 'redmapper')
rmcatfile = os.path.join(rmpath, 'redmapper_isedfit_v5.10_centrals.fits.gz')
rminfo = fitsio.FITS(rmcatfile)
rmcat = rminfo[1].read(columns=['Z', 'MSTAR_50', 'ILUM', 'LAMBDA_CHISQ', 
                                'P_CEN', 'P_SAT', 'MSTAR_AVG'])

def indhist(data, bounds):
    catalog = data
    indices = np.where((catalog >= bounds[0]) * (catalog <= bounds[1]))[0]
    return indices

def getdata(physics):
    for jj in range(numslices):
        physics[jj+1] = indhist(rmcat['MSTAR_AVG'][keep], 
                                [masschunks[1][jj], masschunks[1][jj+1]])
    return;

# Convert dictionary to matrix
def dicttomatrix(dictionary):
    dataarray = []
    for ii in range(numslices):
        dataarray.append(np.array(list(dictionary.items())[ii][1]))
    return dataarray

''' This function accepts two 1x2 arrays. The first is the richness bounds. 
The second is the redshift bounds.'''
def cutter(lbounds, zbounds):
    keep = np.where((rmcat['LAMBDA_CHISQ'] > lbounds[0]) * 
                    (rmcat['LAMBDA_CHISQ'] < lbounds[1]) * 
                    (rmcat['Z'] > zbounds[0]) * 
                    (rmcat['Z'] < zbounds[1]))[0]
    return keep

masses = {}
numslices = 100

keep = cutter([0, 100], [0, 1])
# fig, ax = plt.subplots()
# ax.scatter(rmcat['Z'][keep], rmcat['LAMBDA_CHISQ'][keep])

masssep = np.linspace(min(rmcat['MSTAR_AVG'][keep]), 
                      max(rmcat['MSTAR_AVG'][keep]), numslices)
mdiff = masssep[1]-masssep[0]
# Do I need to index based on mass too?... Make a damn histogram!
masschunks = plt.hist(rmcat['MSTAR_AVG'][keep], bins=numslices)

getdata(masses)
massarray = dicttomatrix(masses)

plt.figure()
for ii in range(numslices):
    if len(massarray[ii] > 0):
        phic = np.sum(rmcat['P_CEN'][massarray[ii]])/(len(massarray[ii])*mdiff)
        phis = np.sum(rmcat['P_SAT'][massarray[ii]])/(len(massarray[ii])*mdiff)
        phi = phic + phis 
        plt.plot(masschunks[1][ii], phi, 'ko')

plt.show()
