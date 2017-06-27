#! usr/bin/env python

import pyfits
import numdisplay

file = '/home/obsastro2/macs1149_f105w.fits.gz'
hdr = pyfits.getheader(file)


im = pyfits.getdata(file,0)
type(im)

im.shape

hdr[0:10]

len(hdr)
 

cutout = im [2000:2600,2000:2500]
pyfits.writeto('junk101.fits',cutout)



