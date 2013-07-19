#! /usr/bin/env python

# montage -label 'MACS0329' macs0329-bcg.png -label 'MACS1206' macs1206-bcg.png junk.png

import sextutils as se
import numpy as np
import pyfits
from astLib import astWCS
from astLib import astImages
from scipy.misc import imread
import glob
import Image
import math

def bcg_cutout(clustername,ra,dec,dirname,width):
    #read image, header
    rootpath = '/moustakas-archive/clash-archive/'+dirname+'/HST/'
    fitspath = rootpath + 'images/mosaicdrizzle_image_pipeline/scale_65mas/'

    pngfile = rootpath + 'color_images/mosaicdrizzle_image_pipeline/'+dirname+'.png'
    fitsfile = glob.glob(fitspath + clustername + '_mosaic_065mas_wfc3ir_f160w_drz_*.fits*')[0]
#    print fitsfile
    hdr = pyfits.getheader(fitsfile)
#   im = pyfits.getdata(fitsfile)
    png = Image.open(pngfile)
#    print png.size
    wcs = astWCS.WCS(hdr,mode='pyfits')
    #slice image

    rawidth = width/math.cos(math.radians(dec))
    decwidth = width
#    print rawidth, decwidth
    ra_left = ra+rawidth
    ra_right = ra-rawidth
    dec_top = dec+decwidth
    dec_bottom = dec-decwidth
#    print ra_left, ra_right, dec_top, dec_bottom

    xy_upperleft = wcs.wcs2pix(ra_left,dec_top)
    xy_bottomright = wcs.wcs2pix(ra_right,dec_bottom)
#    print xy_upperleft, xy_bottomright
#    print xy_upperleft[0]-xy_bottomright[0]
#    print xy_upperleft[1]-xy_bottomright[1]
#    print 2*3600.0*width/0.065
    
#   box = tuple([2400,2300,2700,2600])
    box = tuple([int(xy_upperleft[0]),5000-int(xy_upperleft[1]),
                 int(xy_bottomright[0]),5000-int(xy_bottomright[1])])
#    print box
    cutbcg = png.crop(box)
    
#    cutbcg = astImages.clipImageSectionWCS(im,wcs,ra,dec,width,returnWCS = True)
    return cutbcg
    
if __name__ == "__main__":
    #read BCG file
    bcgfile = '/home/obsastro2/siena-astrophysics/summer2013/clash/clash-sample.txt'    
    bcg = se.se_catalog(bcgfile)    
    outpath = '/home/obsastro2/clash/'
    #loop BCG in BCG file
#    for j in range (1):
#        i = 16
    for i in range (len(bcg.ra)):
        outfile = outpath+str(bcg.shortname[i])+'-bcg.png'
        print i, outfile
                
    #call bcg_cutout
        bcgcut = bcg_cutout(bcg.shortname[i],bcg.ra_bcg[i],bcg.dec_bcg[i],bcg.dirname[i],float(bcg.width[i]/3600.))
        bcgcut.save(outfile, "PNG")
        #write .png
#        mx = float(bcg.stretch[i])
#        astImages.saveBitmap(outfile,cutbcg['data'],[0.0,mx],300,'gray')
