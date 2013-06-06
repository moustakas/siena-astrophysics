#! /usr/bin/env python

import numpy
import numdisplay
import pyfits
from astLib import astWCS
from astLib import astImages

def cutout():
    path = '/home/ioannis/'
    myfile = 'macs1149_f105w.fits.gz'

    im = pyfits.getdata(path+myfile)
    hdr = pyfits.getheader(path+myfile)
    sz = im.shape
    
    wcs = astWCS.WCS(hdr,mode='pyfits')

# get the celestial coordinates at the center
    radec = wcs.pix2wcs(sz[0]/2,sz[1]/2)
    width = 60.0/3600.0 # degrees
    cutim = astImages.clipImageSectionWCS(im,wcs,radec[0],radec[1],width,returnWCS=True)

    astImages.saveBitmap('junk.png',cutim['data'],[0.0,0.5],300,'gray')
