#!/usr/bin/env python

"""
Search for Planet 9 in DECaLS/DR3.

Katie Hoag
2016 June 8
Siena College

"""

import os
import numpy as np
import argparse
import Image, ImageDraw

import pdb  # python debugger

from astropy.io import fits
from astrometry.util.fits import fits_table


def main():

    """ This script creates a set of .jpg files in a candidate_cutouts directory
    which are cutouts of possible Planet Nine candidates obtained from the
    planet9dr3.py script.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--make-cutouts', action='store_true', help='Create cutouts from the DECaLS viewer.')
    parser.add_argument('--make-webpage', action='store_true', help='Create the HTML web content.')
    args = parser.parse_args()

    in_file = ('/home/desi2/candidatesp9/asteroids_decals_dr2.fits')
    out_dir = os.path.join(os.environ.get('HOME'), 'asteroid_cutouts/')

    # --------------------------------------------------------
    # Get thumbnails of objects from DECaLS dr2 and get links to the legacysurvey viewer.
    if args.make_cutouts:
        cand_info = fits_table(in_file)
        # Pre-select asteroids in the ra, dec box you know they exist.
        ramin = 107
        ramax = 130
        decmin = 16
        decmax = 30
        these = np.where((cand_info.ra0>ramin)*(cand_info.ra0<ramax)*
                         (cand_info.dec0>decmin)*(cand_info.dec0<decmax))[0]
        #pdb.set_trace()  # Runs Python Debugger on code up to this line.   
        cand_info = cand_info[these]

        urls = []
        jpgfiles = []
        radius = 15  # radius for putting a circle around central pixel.
        for ii in range(100):
            print('Working on candidate {}'.format(ii))
            ra = cand_info.ra0[ii]
            dec = cand_info.dec0[ii]
        
            jpgurl = 'http://legacysurvey.org/viewer/jpeg-cutout-decals-dr2?ra={:.6f}&dec={:.6f}&pixscale=0.262&size=200'.format(ra, dec)
        
            jpgfile = 'obj-{:03d}.jpg'.format(ii)
            jpgfile = os.path.join(out_dir, jpgfile)
            grab = 'wget --continue -O {:s} "{:s}"' .format(jpgfile, jpgurl)
            print(grab)
            os.system(grab)

            if os.stat(jpgfile).st_size < 18000:  # Remove partial or empty images
                # The cut on filesize takes care of most of the bad images but
                # leaves some behind. If the restriction is any larger,
                # it can remove some valid files.
                os.remove(jpgfile)
            else:
                print(jpgurl)
                jpgfiles.append(jpgfile)
                urls.append(jpgurl)
                # Draw a circle around the object.
                image = Image.open(jpgfile)
                draw = ImageDraw.Draw(image)
                draw.ellipse((100 - radius, 100 - radius, 100 + radius, 100 + radius), outline=66FF00)
                jpg_annot_file = 'obj-{:03d}-annot.jpg'.format(ii)
                jpg_annot_file = os.path.join(out_dir, jpg_annot_file)
                pdb.set_trace()  # Runs Python Debugger on code up to this line.
                
    # ---------------------------------------------------------
    # Create/update the webpage.
    if args.make_webpage:
        
        # for HTML webpage.
        html = open()  #URL?
        html.write('<html><body>/n')
        html.write('<h1> Planet Nine Candidates </h1>/n')
        html.write('<table>/n')
        for thisurl, thisjpg in zip(urls, jpgfiles):
            # add RA, DEC?
            html.write('<tr>/n')
            html.write('<td><a href="{}"><img src="{:s}"></a></td>/n'.format(thisurl, thisjpg))
            html.write('</tr>/n')
        html.write('</table></html>')
        html.close()


if __name__ == '__main__':
    main()
