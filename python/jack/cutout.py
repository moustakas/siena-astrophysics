#! /usr/bin/env python

import numpy as np
import numdisplay
import pyfits
from astLib import astWCS
from astLib import astImages

def cutout():
    path = '/home/obsastro2/'
    myfile = 'macs1149_f105w.fits.gz'

    im = pyfits.getdata(path+myfile)
    hdr = pyfits.getheader(path+myfile)
    sz = im.shape
    
    wcs = astWCS.WCS(hdr,mode='pyfits')

    xcen = np.array([1513.0,2565.0])
    ycen = np.array([2880.0,3076.0])
    width = np.array([30.0,90.0])/3600.0 #degrees

# get the celestial coordinates at the center
    for i in range(len(xcen)):
        radec = wcs.pix2wcs(xcen[i],ycen[i])
        print radec
        ww = float(width[i])
       
        cutim = astImages.clipImageSectionWCS(im,wcs,radec[0],radec[1],ww,returnWCS=True)

        outfile = path+'stamp'+str(i)+'.png'
        print i, outfile
        astImages.saveBitmap(outfile,cutim['data'],[0.0,0.5],300,'gray')

           
