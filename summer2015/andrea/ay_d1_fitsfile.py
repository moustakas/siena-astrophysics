#! /usr/bin/env python

#importing python tools
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt


#setting up the directory path and the file name

dirpath = "/home/work/decam/release/dr1/sweep/"
filename = dirpath + "tractor-sweep-133.fits"


#loads and prints the info and header of the fits file

fits.info(filename)
header_primary =fits.getheader(filename)


#saves the data as a variable and prints the names of the columns

hdu = fits.open(filename)
data = hdu[1].data
print(data.columns)


#attempting to plot a hist the column

NBINS = 100
plt.hist(tbdata['EBV'],NBINS)

