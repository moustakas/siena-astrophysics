#!/usr/bin/env python

'''
Search for Planet 9 in DECaLS/DR3.

Katie Hoag
2016 May 31
Siena College

'''

import os
import numpy as np

import pdb  # python debugger

from astropy.io import fits


def main():

    """ This script creates a set of .jpg files which are cutouts of
    possible Planet Nine candidates obtained from the planet9dr3.py script.
    """

    in_file = os.path.join(os.environ.get('HOME'),
                           'planet9-dr3-candidates.fits')
    out_dir = os.path.join(os.environ.get('HOME'), 'candidate_cutouts/')

    cand_info = fits.getdata(in_file, 1)
    print(cand_info)
    for ii in range(len(in_file)):
        print('Working on candidate {}'.format(ii))
        ra = cand_info['ra']
        dec = cand_info['dec']
        candidate = 
        
        jpeg_url = 'http://legacysurvey.org/viewer/jpeg-cutout-decals-dr2?ra='+ra+'&dec='+dec+'&pixscale=0.100&size=300'
        # fits url?
        # model jpeg url?
        # model fits url?

        os.system('wget "'+jpeg_url+'" -O '+out_dir+ra+dec+'.jpg')
        # os.system for fits
        # os.system for model
        # os.system for fits model

     #pdb.set_trace()  # Runs Python Debugger on code up to this line.


if __name__ == '__main__':
    main()
