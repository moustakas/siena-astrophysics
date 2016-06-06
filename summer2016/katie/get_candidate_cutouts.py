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
from astrometry.util.fits import fits_table


def main():

    """ This script creates a set of .jpg files which are cutouts of
    possible Planet Nine candidates obtained from the planet9dr3.py script.
    """

    in_file = ('/home/desi2/candidatesp9/asteroids_decals_dr2.fits')
    out_dir = os.path.join(os.environ.get('HOME'), 'candidate_cutouts/')

    cand_info = fits_table(in_file)

    urls = []
    for ii in range(20):
        print('Working on candidate {}'.format(ii))
        ra = cand_info.ra0
        dec = cand_info.dec0
        
        jpg_url = 'http://legacysurvey.org/viewer/jpeg-cutout-decals-dr2?ra=%.4f&dec=%.4f&pixscale=0.262&size=200' % (ra, dec)

        urls.append(jpg_url)

        images = []
        for ii,(jpg_url) in enumerate(urls):
            out_jpgs = 'obj-%03i.jpg' % ii
            if not os.path.exists(out_jpgs):
                grab = 'wget --continue -O "%s" "%s"' % (out_jpgs, urls)
                #print(grab)
                os.system(grab)
            images.append(out_jpgs)

    print('<html><body>')
    for ii,(jpg_url,fn) in enumerate(zip(urls,images)):
        print('<a href="http://legacysurvey.org/viewer/?ra=%.4f&dec=%.4f"><img src="%s"></a>' % (ra, dec, fn))
    print('</body></html>')

    pdb.set_trace()  # Runs Python Debugger on code up to this line.


if __name__ == '__main__':
    main()
