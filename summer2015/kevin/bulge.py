#!/usr/bin/env python

import galsim
import os


def bulge():

    bulge_n = 3.5
    bulge_re = 2.3
    disk_n = 1.0         
    disk_r0 = 0.85       
    bulge_frac = 0.3
    
    bulge = galsim.Sersic(bulge_n, half_light_radius=bulge_re)
    disk = galsim.Sersic(disk_n, scale_radius=disk_r0)           
    gal = bulge_frac * bulge + (1-bulge_frac) * disk
    #gal = gal.shear(q=ba, beta=beta*galsim.degrees)
    image = gal.drawImage(scale=0.25)
    image.write('bulge.fits',clobber=True)

if __name__ == "__main__":
    bulge()
