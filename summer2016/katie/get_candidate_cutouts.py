#!/usr/bin/env python

'''
Search for Planet 9 in DECaLS/DR3.

Katie Hoag
2016 May 31
Siena College

'''

import os
import numpy as np
from astropy.io import fits

import pdb


def main():

    ''' This script creates a set of .jpg files which are cutouts of possible Planet Nine candidates obtained from the planet9dr3.py script.
'''

    in_file = os.path.join(os.environ.get('HOME'), 'planet9-dr3-candidates.fits')
    out_dir = os.path.join(in_dir, 'candidate_cutouts/')
    
    cutout_size = 100  # number of pixels per side of the cutout
    
    for ii in range(len(in_file)):
        print('Working on candidate {}'.format(ii))
        ra = 
        dec =
        
        jpeg_url = 'http://legacysurvey.org/viewer/jpeg-cutout-decals-dr2?ra='+ra+'&dec='+dec+'&pixscale=0.262&size=100'
    



if __name__ == '__main__':
    main()
