#!/usr/bin/env python

''' This script will create a set of .jpg files which are cutouts of possible
Planet nine candidates obtained from the planet9dr3.py script

'''

import os
import numpy as np
from astropy.io import fits


def main():

    in_dir = os.path.join(os.environ.get('HOME'), 'planet9-dr3-candidates.fits')
    out_dir = in_dir + 'candidate_cutouts/'
    

    cutout_size = 20  # number of pixels per side of the cutout
    






    
