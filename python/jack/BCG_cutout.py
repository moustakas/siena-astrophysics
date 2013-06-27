#! /usr/bin/env python

import sextutils as se
import numpy as np
import pyfits
from astLib import astWCS
from astLib import astImages
from scipy.misc import imread
#import matplotlib.pyplot as plt
import glob


def bcg_cutout(clustername,ra,dec,dirname,width):
    #read image, header
    rootpath = '/clash-archive/clash_archive/'+dirname+'/HST/'
    fitspath = rootpath + 'images/mosaicdrizzle_image_pipeline/scale_65mas/'

    fitsfile = glob.glob(fitspath + clustername + '_mosaic_065mas_wfc3ir_f160w_drz_*.fits*')[0]
    print fitsfile
    
    im = pyfits.getdata(fitsfile)
    hdr = pyfits.getheader(fitsfile)   
    wcs = astWCS.WCS(hdr,mode='pyfits')
    #width = 30.0/3600.0 #degrees
    #slice image
    cutbcg = astImages.clipImageSectionWCS(im,wcs,ra,dec,width,returnWCS = True)
    return cutbcg
    
if __name__ == "__main__":
    #read BCG file
    bcgfile = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-sample.txt'    
    bcg = se.se_catalog(bcgfile)    
    outpath = '/home/obsastro2/clash/'
    #loop BCG in BCG file
    for i in range (len(bcg.ra)):
        outfile = outpath+str(bcg.shortname[i])+'-bcg.png'
        print i, outfile
                
        #call bcg_cutout
        cutbcg = bcg_cutout(bcg.shortname[i],bcg.ra_bcg[i],bcg.dec_bcg[i],bcg.dirname[i],float(bcg.width[i]/3600.))
        #write .png
        mx = float(bcg.stretch[i])
        astImages.saveBitmap(outfile,cutbcg['data'],[0.0,mx],300,'gray')
