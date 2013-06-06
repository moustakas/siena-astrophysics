#! usr/bin/env python


#file = '/home/obsastro2/macs1149_f105w.fits.gz'
#hdr = pyfits.getheader(file)


#im = pyfits.getdata(file,0)
#type(im)

#im.shape

#hdr[0:10]

#len(hdr)
 

#cutout = im [2000:2600,2000:2500]
#pyfits.writeto('junk101.fits',cutout)

import pyfits

class Snakies(object):
   def  __init__(self,afile,firstx, firsty, lastx, lasty):
        self.afile = afile
        self.firstx = firstx
        self.firsty = firsty
        self.lastx = lastx
        self.lasty = lasty
        

   image = pyfits.getdata(self.afile,0)
   
   def cutout(self):
       cut = self.image[self.firstx,self.lastx,self.firsty,self.lasty]
       pyfits.writeto('crops.fits',cut)
