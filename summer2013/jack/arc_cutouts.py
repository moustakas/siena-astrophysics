#! /usr/bin/env python

import sextutils as sex
import numpy as np
import pyfits
from astLib import astWCS
from astLib import astImages
import argparse
from scipy.misc import imread
import matplotlib.pyplot as plt

def get_cutout(clustername):

    outpath = '/home/obsastro2/clash/'
    datapath = '/home/obsastro2/siena-astrophysics/summer2013/clash/'
    rootpath = '/clash-archive/clash_archive/'+clustername+'/HST/'
    fitspath = rootpath + '/images/mosaicdrizzle_image_pipeline/scale_65mas/'
    colorpath = rootpath + 'color_images/mosaicdrizzle_image_pipeline/'
    
    fitsfile = fitspath + clustername + '_mosaic_065mas_wfc3ir_f160w_drz_20110815.fits.gz'
    colorfile = colorpath + clustername + '_ir.png'

    im = pyfits.getdata(fitsfile)
    hdr = pyfits.getheader(fitsfile)   
    wcs = astWCS.WCS(hdr,mode='pyfits')

    arcfile = datapath+clustername+'-arcs.txt'
    arcs = sex.se_catalog(arcfile)
    
    width = 30.0/3600.0 # degrees
    
    for i in range(len(arcs.ra)):
        cutim = astImages.clipImageSectionWCS(im,wcs,arcs.ra[i],arcs.dec[i],width,returnWCS=True)

        gal = str(arcs.id[i]).split('.')
#       outfile = path+clustername+'-arc'+str(arcs.id[i])+'.png'
        outfile = outpath+clustername+'-arc'+gal[0]+'-im'+gal[1]+'.png'
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
