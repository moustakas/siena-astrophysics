#! /usr/bin/env python

import numpy
import numdisplay
import pyfits
from astLib import astWCS
from astLib import astImages

def cutout():
    path = '/Users/ioannis/archive/macs1206/HST/images/mosaicdrizzle_image_pipeline/scale_65mas/'
    file = 'macs1206_mosaic_065mas_wfc3ir_f105w_drz_20110815.fits.gz'

    im = pyfits.getdata(path+file)
    hdr = pyfits.getheader(path+file)
    sz = im.shape
    
    wcs = astWCS.WCS(hdr,mode='pyfits')

# get the celestial coordinates at the center
    radec = wcs.pix2wcs(sz[0]/2,sz[1]/2)
    width = 20.0/3600.0 # degrees
    cutim = astImages.clipImageSectionWCS(im,wcs,radec[0],\
                                              radec[1],width,returnWCS=True)
    out = numpy.array(cutim['clippedSection'])

    astImages.saveBitmap('junk.png',out,'smart','1024','grey')

    astImages.saveBitmap('junk.fits',im,'smart','1024','grey')
    

    
    cutout = astImages.clipImageSectionPix(im,2500,2500,[200,200])
    pyfits.writeto('junk.fits',fim,hdr) # create new fits file

    saveFITS('junk.fits',cutim,imageWCS=wcs)
