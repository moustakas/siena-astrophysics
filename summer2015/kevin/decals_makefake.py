#!/usr/bin/env python

"""Insert artificial (fake) galaxies into the DECaLS imaging and reprocess
through The Tractor.

TODO:
* Model a mixture of object types.
"""
from __future__ import division, print_function

import os
import argparse

import galsim
import numpy as np

from tractor import psfex
from tractor.basics import GaussianMixtureEllipsePSF

# Global variables.
decals_dir = os.getenv('DECALS_DIR')
fake_decals_dir = os.getenv('FAKE_DECALS_DIR')

def get_brickinfo(brickname=None):
    """Get info on this brick and on the CCDs touching it.

    """
    from astropy.io import fits
    from projects.desi.common import *

    allbrickinfo = fits.getdata(os.path.join(decals_dir,'decals-bricks.fits'),1)
    allccdinfo =  fits.getdata(os.path.join(decals_dir,'decals-ccds-zeropoints.fits'))

    brickinfo = allbrickinfo[np.where((allbrickinfo['brickname']==brickname)*1==1)]

    # Get all the CCDs that touch this brick
    decals = Decals()
    wcs = wcs_for_brick(decals.get_brick_by_name(brickname))
    these = ccds_touching_wcs(wcs,decals.get_ccds())
    #ccdinfo = decals.ccds_touching_wcs(targetwcs)
    ccdinfo = allccdinfo[these]

    return brickinfo, wcs, ccdinfo


def build_priors(nobj=20,brickname=None,objtype='ELG',ra_range=None,dec_range=None):
    """Choose priors according to the type of object.  Will eventually generalize
       this so that a mixture of object types can be simulated.

    """

    # Assign central coordinates uniformly
    ra = np.random.uniform(ra_range[0],ra_range[1],nobj)
    dec = np.random.uniform(dec_range[0],dec_range[1],nobj)

    if objtype.upper()=='ELG':
        # Disk parameters
        disk_n_range = [1.0,1.0]
        disk_r50_range = [0.5,2.5]
        disk_ba_range = [0.2,1.0]

        disk_n = np.random.uniform(disk_n_range[0],disk_n_range[1],nobj)
        disk_r50 = np.random.uniform(disk_r50_range[0],disk_r50_range[1],nobj)
        disk_ba = np.random.uniform(disk_ba_range[0],disk_ba_range[1],nobj)
        disk_phi = np.random.uniform(0,180,nobj)

        ## Bulge parameters
        #bulge_r50_range = [0.1,1.0]
        #bulge_n_range = [3.0,5.0]
        #bdratio_range = [0.0,1.0] # bulge-to-disk ratio

        # Magnitudes and colors
        rmag_range = [17.0,19.0]
        gr_range = [-0.3,0.5]
        rz_range = [0.0,1.5]

        rmag = np.random.uniform(rmag_range[0],rmag_range[1],nobj)
        gr = np.random.uniform(gr_range[0],gr_range[1],nobj)
        rz = np.random.uniform(rz_range[0],rz_range[1],nobj)

        # Pack into a Table.

        priors = fits.BinTableHDU.from_columns([
            fits.Column(name='ra',format='f4',array=ra),
            fits.Column(name='dec',format='f4',array=dec),
            fits.Column(name='r',format='f4',array=rmag),
            fits.Column(name='gr',format='f4',array=gr),
            fits.Column(name='rz',format='f4',array=rz),
            fits.Column(name='disk_n',format='f4',array=disk_n),
            fits.Column(name='disk_r50',format='f4',array=disk_r50),
            fits.Column(name='disk_ba',format='f4',array=disk_ba),
            fits.Column(name='disk_phi',format='f4',array=disk_phi)])
        priors.writeto(fake_decals_dir+brickname+'_priors.fits',clobber=True)


def copyfiles():
    """Copy the CP-processed images, inverse variance maps, and bad-pixel masks we
    need from DECALS_DIR to FAKE_DECALS_DIR, creating directories as necessary.

    """
    import shutil
    for cpimage in list(set(allccds.cpimage)):
        outcpdir = os.path.join(fake_decals_dir,'images','decam',cpimage.split('/')[1])
        if not os.path.isdir(outcpdir):
            os.makedirs(outcpdir)
            outfile = os.path.join(outcpdir,cpimage.split('/')[2]).strip()
            print('Copying image {}'.format(outfile))
            shutil.copyfile(os.path.join(decals_dir,'images',cpimage).strip(),outfile)

            outfile = os.path.join(outcpdir,cpimage.split('/')[2]).replace('ooi','oow').strip()
            print('Copying weight map {}'.format(outfile))
            shutil.copyfile(os.path.join(decals_dir,'images',
                                         cpimage).strip().replace('ooi','oow'),outfile)
                                         
            outfile = os.path.join(outcpdir,cpimage.split('/')[2]).replace('ooi','ood').strip()
            print('Copying bad pixel mask {}'.format(outfile))
            shutil.copyfile(os.path.join(decals_dir,'images',
                                         cpimage).strip().replace('ooi','ood'),outfile)
                


def insert_simobj():

    gsparams = galsim.GSParams(maximum_fft_size=2L**30L)

    stampwidth = 95 # postage stamp width [pixels, roughly 14 arcsec]
    stampbounds = galsim.BoundsD(0,stampwidth,0,stampwidth)
    
    band = np.array(['g','r','z'])

    # Calculate the g, r, and z band fluxes and stack them in an array.
    gflux = 10**(-0.4*(rmag+grmag)) # [nanomaggies]
    rflux = 10**(-0.4*rmag)
    zflux = 10**(-0.4*(rmag-rzmag))
    flux = np.vstack([gflux,rflux,zflux])
        
    for ccd in allccds:
        print(ccd.ra,ccd.dec)

        # Get the filenames we need.
        cpname = ccd.cpimage.strip()
        calname = ccd.calname.strip()
        thisband = np.where((ccd.filter==band)*1)
        
        imfile = os.path.join(fake_decals_dir,'images',cpname)
        ivarfile = imfile.replace('ooi','oow')
        psffile = os.path.join(decals_dir,'calib','decam','psfex',calname+'.fits')
        wcsfile = os.path.join(decals_dir,'calib','decam','astrom-pv',calname+'.wcs.fits')
        
        # Read the data.
        print('Reading extension {} of image {}'.format(ccd.ccdnum,imfile))
        im = galsim.fits.read(imfile,hdu=ccd.ccdnum)       # [ADU]
        invvar = galsim.fits.read(ivarfile,hdu=ccd.ccdnum) # [1/ADU^2]
        hdr = galsim.fits.FitsHeader(imfile,hdu=ccd.ccdnum)
        ivarhdr = galsim.fits.FitsHeader(ivarfile,hdu=ccd.ccdnum)
        
        # Get the WCS info and initialize the PSF
        wcs, origin = galsim.wcs.readFromFitsHeader(galsim.fits.FitsHeader(wcsfile))
        initpsf = psfex.PsfEx(psffile,im.xmax,im.ymax,ny=13,nx=7,
                              psfClass=GaussianMixtureEllipsePSF,K=2)

        # Get the zeropoint info for this image and CCD
        zptinfo = zpts[np.where(((zpts['CPIMAGE']==ccd.cpimage)*1)*((zpts['EXTNAME']==ccd.extname)*1))]
        magzpt = float(zptinfo['CCDZPT'] + 2.5*np.log10(zptinfo['EXPTIME']))
        gain = float(zptinfo['ARAWGAIN']) # [electron/ADU]

        # Loop on each fake galaxy.
        for iobj in range(nobj):
            # Get the position of the galaxy on the CCD and build the PSF.
            pos = wcs.toImage(galsim.CelestialCoord(
                ra[iobj]*galsim.degrees,dec[iobj]*galsim.degrees))
            xpos = int(pos.x)
            ypos = int(pos.y)
            print(iobj, xpos, ypos, ndisk[iobj], disk_r50[iobj])

            galflux = float(flux[thisband,iobj]*10**(0.4*magzpt)) # [ADU]

            wcslocal = wcs.local(image_pos=pos)
            offset = galsim.PositionD(pos.x-xpos,pos.y-ypos)
            pixscale, shear, theta, flip = wcslocal.getDecomposition() # get the pixel scale

            psfim = PsfEx.instantiateAt(initpsf,xpos,ypos)[5:-5,5:-5]
            psf = galsim.InterpolatedImage(galsim.Image(psfim),scale=pixscale,flux=1.0)
            psf_centroid = psf.centroid()
            psf = psf.shift(-psf_centroid.x,-psf_centroid.y)

            # Build the postage stamp of the object convolved with the PSF.
            #bulge = galsim.Sersic(n=nbulge[iobj],half_light_radius=ndisk[iobj],
            #                      gsparams=gsparams,flux=flux[iband,iobj])
            #disk = galsim.Sersic(ndisk[iobj],scale_radius=disk_r50[iobj])
            #stamp = bulge_frac[iobj] * bulge + (1-bulge_frac[iobj]) * disk
            gal = galsim.Sersic(ndisk[iobj],scale_radius=disk_r50[iobj],flux=galflux,gsparams=gsparams)
            #gal = galsim.Sersic(ndisk[iobj],scale_radius=disk_r50[iobj],
            #                    flux=float(flux[thisband,iobj]*10**(-0.4*magzpt)))
            gal = gal.shear(q=axisratio[iobj],beta=phi[iobj]*galsim.degrees)
            gal = galsim.Convolve([gal,psf])
            #gal = gal.shift(pos.x-xpos,pos.y-ypos) # apply sub-pixel shift
            
            #stamp = gal.drawImage(wcs=wcslocal,method='no_pixel')
            stamp = gal.drawImage(nx=stampwidth,ny=stampwidth,wcs=wcslocal,method='no_pixel')
            stamp.setCenter(xpos,ypos)
            
            overlap = stamp.bounds & im.bounds
            if (overlap.xmax>=0 and overlap.ymax>=0 and overlap.xmin<=im.xmax and
                overlap.ymin<=im.ymax and overlap.area()>0):
                print('Adding object {} to x={}, y={}, ra={}, dec={} with flux {}'.format(
                        iobj,pos.x,pos.y,ra[iobj],dec[iobj],galflux))
                    
                stamp = stamp[overlap]

                # Add Poisson noise
                varstamp = invvar[overlap].copy()
                varstamp.invertSelf() # [ADU^2]
                medvar = np.median(varstamp.array[varstamp.array>0])
                neg = np.where(varstamp.array<(0.2*medvar))
                if neg[0].size>0:
                    varstamp.array[neg] = medvar
                        
                stamp *= gain # [electrons]
                varstamp *= gain**2
                varstamp += stamp
                
                stamp.addNoise(galsim.VariableGaussianNoise(galsim.BaseDeviate(),varstamp))

                stamp /= gain         # [ADU]
                varstamp /= gain**2   # [ADU^2]
                varstamp.invertSelf() # [1/ADU^2]

                im[overlap] += stamp
                invvar[overlap] += varstamp

            # Writes the images to the output directory.
            outfile = os.path.join(fake_decals_dir,'images',ccd.cpimage.strip())
            print('Updating extension {} of image {}'.format(ccd.ccdnum,outfile))
            fits.update(imfile,im.array,ext=ccd.ccdnum,header=fits.Header(hdr.items()))
            fits.update(ivarfile,invvar.array,ext=ccd.ccdnum,header=fits.Header(ivarhdr.items()))

            #galsim.fits.write(im,file_name='junk.fits',clobber=True)

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='DECaLS simulations.')
    parser.add_argument('--nobj', type=long, default=None, metavar='', 
                        help='number of objects to simulate (required input)')
    parser.add_argument('--objtype', type=str, default='ELG', metavar='', 
                        help='object type (ELG, LRG, PSF, MIX)') 
    parser.add_argument('--brick', type=str, default='2428p117', metavar='', 
                        help='simulate objects in this brick')
    parser.add_argument('--zoom', nargs=4, type=int, metavar='', 
                        help='see runbrick.py (default is to populate the full brick)')
    #parser.add_argument('--build-priors', action='store_false', 
    #                    help="""Build the object priors.""")

    args = parser.parse_args()
    if args.nobj is None:
        parser.print_help()
        sys.exit(1)

    objtype = args.objtype.upper()
    brickname args.brick

    print('Working on brick {}'.format(brickname))
        
    # Get the brick and CCD info
    brickinfo, brickwcs, ccdinfo = get_brickinfo(brickname)

    if nargs.zoom is None:
        ra_range = [info['ra1'],info['ra2']]
        dec_range = [info['dec1'],info['dec2']]
    else:
        pixscale = 0.262/3600.0 # average pixel scale [deg/pixel]
        zoom = nargs.zoom
        dx = zoom[1]-zoom[0]
        dy = zoom[3]-zoom[2]

        ra, dec = wcs.pixelxy2radec(zoom[0]+dx/2,zoom[2]+dy/2)
        ra_range = [ra-dx*pixscale/2,ra+dx*pixscale/2]
        dec_range = [ra-dy*pixscale/2,dec+dy*pixscale/2]
    
    # Build the prior parameters.
    build_priors(args.nobj,brickname,objtype,ra_range,dec_range)

    # Simulate each individual object.
    insert_simobj()


if __name__ == "__main__":
    main()
