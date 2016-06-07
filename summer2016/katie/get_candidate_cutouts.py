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
    out_dir = os.path.join(os.environ.get('HOME'), 'asteroid_cutouts/')

    cand_info = fits_table(in_file)
    # Pre-select asteroids in the ra, dec box you know they exist.
    ramin = 107
    ramax = 130
    decmin = 16
    decmax = 30
    these = np.where((cand_info.ra0>ramin)*(cand_info.ra0<ramax)*
                     (cand_info.dec0>decmin)*(cand_info.dec0<decmax))[0]
    cand_info = cand_info[these]

    urls = []
    jpgfiles = []
    for ii in range(4):
        print('Working on candidate {}'.format(ii))
        ra = cand_info.ra0[ii]
        dec = cand_info.dec0[ii]
        
        jpgurl = 'http://legacysurvey.org/viewer/jpeg-cutout-decals-dr2?ra={:.4f}&dec={:.4f}&pixscale=0.262&size=200'.format(ra, dec)
        
        jpgfile = 'obj-{:03d}.jpg'.format(ii)
        jpgfile = os.path.join(out_dir, jpgfile)
        grab = 'wget --continue -O {:s} {:s}' .format(jpgfile, jpgurl)
        print(grab)
   
        os.system(grab)
        if os.stat(jpgfile).st_size == 0:
            os.remove(jpgfile)
        else:
            print(jpgurl)
            jpgfiles.append(jpgfile)
            urls.append(jpgurl)
    pdb.set_trace()  # Runs Python Debugger on code up to this line.
    print('<html>')
    print('<head> Planet Nine Candidates </head>')
    print('<body>')
    for thisurl, thisjpg in zip(urls, jpgfiles):
        print('<div class="image">')
        print('<a href="{}"><img src="{:s}"></a>'.format(thisurl, thisjpg))
        print('<div class="caption"> Image of {:s} </div>' .format(thisjpg))
        print('</div>')
    print('</body></html>')



if __name__ == '__main__':
    main()
