#!/usr/bin/env python

"""
Search for Planet 9 in DECaLS/DR2.

Katie Hoag
2016 June 9
Siena College

"""

import os
import numpy as np
import argparse
import Image, ImageDraw

import pdb

from astropy.io import fits
from astrometry.util.fits import fits_table


def main():

    """ This script creates a set of .jpg files in a candidate_cutouts directory
    which are cutouts of possible Planet 9 candidates obtained from the DECaLS
    DR2 Tractor files.

    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--make-cutouts', action='store_true', help='Create cutouts from the DECaLS viewer.')
    parser.add_argument('--make-webpage', action='store_true', help='Create or remake the HTML web content.')
    args = parser.parse_args()

    in_file = os.path.join(os.environ.get('HOME'), 'candidatesp9/planet9-dr2-candidates.fits')
    out_dir = os.path.join(os.environ.get('HOME'), 'candidatesp9/candidate_cutouts/')

    # Read the sample.
    cand_info = fits_table(in_file)
    #pdb.set_trace()  # Runs Python Debugger on code up to this line.   

    # --------------------------------------------------------
    # Get thumbnails of objects from DECaLS dr2 and get links to the legacysurvey viewer.
    if args.make_cutouts:
        urls = []
        jpgfiles = []
        jpgannotfiles = []        

        radius = 15  # radius for putting a circle around central pixel.
        for ii in range(len(cand_info)):
            print('Working on candidate {}'.format(ii))
            ra = cand_info.ra[ii]
            dec = cand_info.dec[ii]
        
            jpgurl = 'http://legacysurvey.org/viewer/jpeg-cutout-decals-dr2?ra={:.6f}&dec={:.6f}&pixscale=0.262&size=200'.format(ra, dec)
            jpgfile = 'cand-{:04d}.jpg'.format(ii)
            jpgfile = os.path.join(out_dir, jpgfile)
            grab = 'wget --continue -O {:s} "{:s}"' .format(jpgfile, jpgurl)
            os.system(grab)

            if os.stat(jpgfile).st_size < 0:  # Remove empty images
                os.remove(jpgfile)
            else:
                jpgfiles.append(jpgfile)
                urls.append(jpgurl)
                # Draw a circle around the object.
                image = Image.open(jpgfile)
                draw = ImageDraw.Draw(image)
                draw.ellipse((100 - radius, 100 - radius, 100 + radius, 100 + radius), outline="#66FF00")
                jpg_annot_file = image.save('/home/desi2/candidatesp9/candidate_cutouts/cand-{:04d}-annot.jpg'.format(ii), 'JPEG')
                jpgannotfiles.append(jpg_annot_file)
                #pdb.set_trace()  # Runs Python Debugger on code up to this line.
                
    # ---------------------------------------------------------
    # Create/update the webpage.
    if args.make_webpage:
        
        # for HTML webpage.
        html = open("/home/desi2/candidatesp9/index.html", 'w')
        html.write('<html><body>\n')
        html.write('<a name="top"></a>\n')
        html.write('<h1> Possible Planet 9 Candidates from DR2 with Z </h1>\n')
        html.write('<h4> These candidates were obtained by passing data through strict requirements. Candidates are not yet matched against known asteroids and objects.</h4>')
        html.write('<table border="1" style="width:60%">\n')

        for ii in range(len(cand_info)):
        #for ii in range(100):
            ra = cand_info.ra[ii]
            dec = cand_info.dec[ii]
        
            jpgurl = 'http://legacysurvey.org/viewer/jpeg-cutout-decals-dr2?ra={:.6f}&dec={:.6f}&pixscale=0.262&size=200'.format(ra, dec)
            viewerurl = 'http://legacysurvey.org/viewer/?ra={:.6f}&dec={:.6f}&zoom=16'.format(ra,dec)
        
            jpgfile = os.path.join('candidate_cutouts/', 'cand-{:04d}-annot.jpg'.format(ii))
            if os.path.exists('/home/desi2/candidatesp9/'+jpgfile):
                html.write('<tr>\n')
                html.write('<td align=center>Image {:d}<br/>Candidate with RA={:.6f} <br/>and DEC={:.6f} </td>\n'.format(ii+1, ra, dec))
                html.write('<td align=center><a href="{}"><img src="{:s}"></a></td>\n'.format(jpgurl, jpgfile))
                html.write('<td align=center><a href="{}"> LegacySurvey Viewer </a></td>\n'.format(viewerurl))
                html.write('</tr>\n')
            if ii % 10 == 0 and ii > 0:
                html.write('<td colspan="3" align=right><a href="#top"> <font size=4> Back to Top of Page </font> </a></td>\n')
        html.write('</table>\n')
        html.write('</html>\n')
        html.close()


if __name__ == '__main__':
    main()
