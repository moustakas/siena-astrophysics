#! /usr/bin/env python

import pyfits
import numdisplay

pyfits.open("macs1149_f105w.fits.gz")

header_primary = pyfits.getheader("macs1149_f105w.fits.gz")

file = '/home/obsastro2/macs1149_f105w.fits.gz'

type(file)
file.shape


