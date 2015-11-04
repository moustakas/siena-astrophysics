#!/usr/bin/env python

# Build the multi-color catalogs for Ben's HFF project

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
#from astrometry.libkd import spherematch_c

def build_cat(cluster,topdir):
    
    cat = ascii.read(topdir+cluster+'_multicolor_nir.cat',
                     format='sextractor')

    f105w = cat['hst_wfc3_ir_f105w_mag_bpz']
    f140w = cat['hst_wfc3_ir_f140w_mag_bpz']
    f160w = cat['hst_wfc3_ir_f160w_mag_bpz']

    ascii.write(cat,topdir+cluster+'_hff.cat',
                include_names=['number',
                               'alpha_j2000',
                               'delta_j2000',
                               'hst_acs_wfc_f814w_mag_bpz',
                               'hst_wfc3_ir_f105w_mag_bpz',
                               'hst_wfc3_ir_f140w_mag_bpz',
                               'hst_wfc3_ir_f160w_mag_bpz'
                           ],
                names=['number','ra','dec','f814w','f105w','f140w','f160w'],
                dtype=['i8','f8','f8','f8','f8','f8','f8'],
                format='basic')
    
if __name__ == '__main__':

    topdir = '/Users/ioannis/research/people/ben/'
    for cat in ['a2744','m0416']:
        build_cat(cat,topdir)
