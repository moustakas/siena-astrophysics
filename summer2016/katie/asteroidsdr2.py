#!/usr/bin/env python

"""
Search for Asteroids in DECaLS/DR2.

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

    """ This script creates a set of .jpg files in an asteroid_cutouts directory
    which are cutouts of asteroids obtained from the asteroids_decals_dr2/fits
    file.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--make-cutouts', action='store_true', help='Create cutouts from the DECaLS viewer.')
    parser.add_argument('--make-webpage', action='store_true', help='Create the HTML web content.')
    args = parser.parse_args()

    in_file = ('/home/desi2/candidatesp9/asteroids_decals_dr2.fits')
    out_dir = os.path.join(os.environ.get('HOME'), 'asteroids/asteroid_cutouts/')
  

    # Read the sample and cull it.
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

    # --------------------------------------------------------
    # Get thumbnails of objects from DECaLS dr2 and get links to the legacysurvey viewer.
    if args.make_cutouts:
        urls = []
        jpgfiles = []
        jpgannotfiles = []        

        radius = 15  # radius for putting a circle around central pixel.
        for ii in range(len(cand_info)):
            print('Working on candidate {}'.format(ii))
            ra = cand_info.ra0[ii]
            dec = cand_info.dec0[ii]
        
            jpgurl = 'http://legacysurvey.org/viewer/jpeg-cutout-decals-dr2?ra={:.6f}&dec={:.6f}&pixscale=0.262&size=200'.format(ra, dec)
            jpgfile = 'obj-{:04d}.jpg'.format(ii)
            jpgfile = os.path.join(out_dir, jpgfile)
            grab = 'wget --continue -O {:s} "{:s}"' .format(jpgfile, jpgurl)
            os.system(grab)

            if os.stat(jpgfile).st_size < 18000:  # Remove partial or empty images
                # The cut on filesize takes care of most of the bad images but
                # leaves some behind. If the restriction is any larger,
                # it can remove some valid files.
                os.remove(jpgfile)
            else:
                jpgfiles.append(jpgfile)
                urls.append(jpgurl)
                # Draw a circle around the object.
                image = Image.open(jpgfile)
                draw = ImageDraw.Draw(image)
                draw.ellipse((100 - radius, 100 - radius, 100 + radius, 100 + radius), outline="#66FF00")
                jpg_annot_file = image.save('/home/desi2/asteroids/asteroid_cutouts/obj-{:04d}-annot.jpg'.format(ii), 'JPEG')
                jpgannotfiles.append(jpg_annot_file)
                #pdb.set_trace()  # Runs Python Debugger on code up to this line.
                
    # ---------------------------------------------------------
    # Create/update the webpage.
    if args.make_webpage:
        
        # for HTML webpage.
        html = open("/home/desi2/asteroids/index.html", 'w')
        html.write('<html><body>\n')
        html.write('<a name="top"></a>\n')
        html.write('<h1> Asteroids in DR2 </h1>\n')
        html.write('<h3> Location of Asteroids is within a range of RA 107:130 degrees and DEC 16:30 degrees </h3>')
        html.write('<table border="1" style="width:50%">\n')

        for ii in range(len(cand_info)):
        #for ii in range(100):
            ra = cand_info.ra0[ii]
            dec = cand_info.dec0[ii]
        
            jpgurl = 'http://legacysurvey.org/viewer/jpeg-cutout-decals-dr2?ra={:.6f}&dec={:.6f}&pixscale=0.262&size=200'.format(ra, dec)
            viewerurl = 'http://legacysurvey.org/viewer/?ra={:.6f}&dec={:.6f}&zoom=16'.format(ra,dec)
        
            jpgfile = os.path.join('asteroid_cutouts', 'obj-{:03d}-annot.jpg'.format(ii))
            if os.path.exists('/home/desi2/asteroids/'+jpgfile):
                html.write('<tr>\n')
                html.write('<td align=center>Image {:d}<br/>Asteroid with RA={:.6f} <br/>and DEC={:.6f} </td>\n'.format(ii+1, ra, dec))
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
