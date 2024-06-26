#! /usr/bin/env python

from astropy.io import ascii
from astropy.io import fits
from astropy import wcs
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    filt = ['f814w','f125w','f140w','f160w']
    #filt = ['f435w','f606w','f814w','f105w','f125w','f140w','f160w']
    nfilt = len(filt)

    width = 55 # pixels
    non_linear = 0.1

    # Read the master catalog
    codedir = '/Users/ioannis/repos/git/siena-astrophysics/research/hff/'
    cat = ascii.read(codedir+'z9_candidates.cat',format='sextractor')
    cluster = np.unique(cat['cluster'])

    topdir = '/Users/ioannis/research/people/ben/'

    for cl in cluster:
        allobj = cat[np.where((cl==cat['cluster'])*1)]
        for obj in allobj:
            fig = plt.figure(figsize=(3*nfilt,3))
            plt.subplots_adjust(wspace=0.02,bottom=0.02,top=1,left=0,right=1)
            #fig, ax = plt.subplots(nrows=1,ncols=nfilt,sharey=True,sharex=True)
            #plt.setp([a.get_yticklabels() for a in fig.axes],visible=False)
            #plt.setp([a.get_xticklabels() for a in fig.axes],visible=False)
            nf = 1
            for ff in filt:
                fitsfile = glob(topdir+'*'+cluster2+'*'+ff+'*')
                print(fitsfile)
                hdulist = fits.open(fitsfile[0])
                # Parse the WCS keywords in the primary HDU
                w = wcs.WCS(hdulist[0].header)
                # w.wcs.print_contents()
                ra = cat[ii]['RA']
                dec = cat[ii]['DEC']
                xycen = w.wcs_world2pix(ra,dec,1)
                xcen = np.round(xycen[0])
                ycen = np.round(xycen[1])
                print(ra, dec, xcen, ycen)
            
            im = hdulist[0].data
            cutout = im[ycen-width/2:ycen+width/2,xcen-width/2:xcen+width/2]

            # get the arcsinh image scaling
            scale_min = np.array(cutout.min())
            #scale_min[scale_min<0] = 0
            scale_max = cutout.max()
            factor = np.arcsinh((scale_max - scale_min)/non_linear)
            #print scale_min, scale_max, factor
            cutout_scaled = np.arcsinh((cutout-scale_min)/non_linear)/factor
            #print cutout_scaled.min(), cutout_scaled.max()
            
            #cutout_scaled = cutout

            plt.subplot(1,nfilt,nf)
            plt.imshow(cutout_scaled,vmin=0.0,vmax=1.0,cmap=plt.get_cmap('gray_r'))
            plt.axis('off')
            plt.text(0.1,0.85,ff.upper(),transform=plt.gca().transAxes,
                     horizontalalignment='left',color='black',fontweight='bold',
                     fontsize=18)
            circ = plt.Circle((width/2,width/2),radius=5,color='g',fill=False,lw=3)
            plt.gca().add_patch(circ)

            if nf==nfilt:
                plt.text(0.9,0.1,cat[ii]['ID'],transform=plt.gca().transAxes,
                         horizontalalignment='right',color='black',fontweight='bold',
                         fontsize=20)
                nf += 1

        pngfile = topdir+'cutouts/'+cluster+'-'+cat[ii]['ID']+'-cutout.png'
        print pngfile
        plt.savefig(pngfile)
        plt.close(fig)
                
        hdulist.close()

