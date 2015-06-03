#!/usr/bin/env python

import galsim
import os


def bulge():

    bulge = galsim.Sersic(bulge_n, half_light_radius=bulge_re)
    disk = galsim.Sersic(disk_n, scale_radius=disk_r0)           
    gal = bulge_frac * bulge + (1-bulge_frac) * disk
    gal = gal.shear(q=ba, beta=beta*galsim.degrees)

if __name__ == "__main__":
    bulge()
