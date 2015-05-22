#! /usr/bin/env python

#importing python tools
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

#####
#importing the sdss files
from sklearn.neighbors import KNeighborsRegressor
from astroML.plotting import scatter_contour


#setting up the directory path and the file name

dirpath = "/home/desi2/"
filename = dirpath + "stars_qsos_sdss.fits"


#loads and prints the info and header of the fits file

fits.info(filename)
header_primary =fits.getheader(filename)


#saves the data as a variable and prints the names of the columns

hdu = fits.open(filename)
data = hdu[1].data
print(data.columns)

#####
#uses the fits file for the data
#data = data[::10]  # truncate for plotting

# Extract colors and spectral class
ug = data['u'] - data['g']
gr = data['g'] - data['r']
# class = data['class']


star = (data['class'] == 'STAR')
qsos = (data['class'] == 'QSO')


#------------------------------------------------------------
# Prepare plot
fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlim(-0.5, 2.5)
ax.set_ylim(-0.5, 1.5)

ax.plot(ug[star], gr[star], '.', ms=4, c='b', label='stars')
ax.plot(ug[qsos], gr[qsos], '.', ms=4, c='r', label='qsos')

ax.legend(loc=2)

ax.set_xlabel('$u-g$')
ax.set_ylabel('$g-r$')

print('Writing test_sdss.png')
plt.savefig('/home/desi2/test_sdss.png')

#!display test_sdss.png &




