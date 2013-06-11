#! /usr/bin/env python

import sextutils as sex
import numpy as np
import pyfits
from astLib import astWCS
from astLib import astImages
import argparse

def get_cutout(clustername):
    
    path = '/Users/ioannis/tmp/'
    imfile = clustername+'_f160w.fits.gz'

    im = pyfits.getdata(path+imfile)
    hdr = pyfits.getheader(path+imfile)
    sz = im.shape    
    wcs = astWCS.WCS(hdr,mode='pyfits')

    arcfile = path+clustername+'-arcs.sex'
    arcs = sex.se_catalog(arcfile)
    
#   xcen = np.array([1513.0,2565.0])
#   ycen = np.array([2880.0,3076.0])
#   width = np.array([30.0,90.0])/3600.0 #degrees
    width = 30.0/3600.0 # degrees
    
# get the celestial coordinates at the center
    for i in range(10):
#   for i in range(len(arcs.ra)):
#      radec = wcs.pix2wcs(xcen[i],ycen[i])
#       ww = float(width[i])
       
        cutim = astImages.clipImageSectionWCS(im,wcs,arcs.ra[i],arcs.dec[i],width,returnWCS=True)

        gal = str(arcs.id[i]).split('.')
#       outfile = path+clustername+'-arc'+str(arcs.id[i])+'.png'
        outfile = path+clustername+'-arc'+gal[0]+'-im'+gal[1]+'.png'
        print i, outfile
        astImages.saveBitmap(outfile,cutim['data'],[0.0,0.5],300,'gray')

           
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Retrieve cluster cutouts.')
    parser.add_argument('cluster', type=str, default=None, help='Cluster name')

    args = parser.parse_args()

    if args.cluster is None:
        print "Need a cluster name!"
        parser.print_help()

    clustername = args.cluster
    get_cutout(clustername)
