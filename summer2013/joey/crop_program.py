#! usr/bin/env python

import pyfits
import numpy
from astLib import astWCS
from astLib import astImages

def crop_program():
   path = '/home/obsastro1/'
   myfile = 'macs1149_f105 w.fits.gz'

   im = pyfits.getdata(path+myfile)


#class Snakies(object):
#   def  __init__(self,afile,firstx, firsty, lastx, lasty):
#        self.afile = afile
#        self.firstx = firstx
#        self.firsty = firsty
#        self.lastx = lastx
#        self.lasty = lasty
#        
#
#   image = pyfits.getdata(self.afile,0)
#   
#   def cutout(self):
#       cut = self.image[self.firstx,self.lastx,self.firsty,self.lasty]
#       pyfits.writeto('crops.fits',cut)
