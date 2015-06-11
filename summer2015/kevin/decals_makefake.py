#!/usr/bin/env python

import galsim
import numpy as np
import os
import shutil
import math
from astropy.io import fits
from projects.desi.common import *
from tractor import psfex
from tractor.basics import GaussianMixtureEllipsePSF

def decals_makefake():

    gsparams = galsim.GSParams(maximum_fft_size=2**16)

    # Set parameter ranges.
    nfake = 1
    min_disk_r50 = 0.5
    max_disk_r50 = 3.0
    min_bulge_r50 = 0.1
    max_bulge_r50 = 1.0
    nbulgemin = 3.0
    nbulgemax = 5.0
    ndiskmin = 0.8
    ndiskmax = 2.0
    min_bulge_frac = 0.0
    max_bulge_frac = 1.0
    phimin = 0.0
    phimax = 180.0
    axismin = 0.2
    axismax = 1.0
    rmagmin = 16.0
    rmagmax = 18.0
    grmin = -0.3
    grmax = 0.5
    rzmin = 0.0
    rzmax = 1.5
   
    imw = 2046
    imh = 4094
    
    # Set directories.
    decals_dir = os.getenv('DECALS_DIR')
    fake_decals_dir = os.getenv('FAKE_DECALS_DIR')
    brickinfo = fits.getdata(decals_dir+'/decals-bricks.fits',1)
    ccdinfo = fits.getdata(decals_dir+'/decals-ccds.fits',1)
    
    # Create an array for both the bricks and bands.
    mybricks = ['2428p117']
    #mybricks = ['2426p197','2428p117']
    band = ['g','r','z']

    # initialize some classes we'll need
    decals = Decals()
    
    #  Loop randomly generating priors.
    for thisbrick in mybricks:
        print('Getting info on brick {}'.format(thisbrick))
        info = brickinfo[np.where((brickinfo['brickname']==thisbrick)*1==1)]
        
        # Randomly generate location and flux parameters.
        ra = np.random.uniform(info['ra1'],info['ra2'],nfake)
        dec = np.random.uniform(info['dec1'],info['dec2'],nfake)
        rmag = np.random.uniform(rmagmin,rmagmax,nfake)
        grmag = np.random.uniform(grmin,grmax,nfake)
        rzmag = np.random.uniform(rzmin,rzmax,nfake)

        # Calculate the g, r, and z band fluxes and stack them in an array.
        gflux = 10**(-0.4*((rmag+grmag)-22.5))
        rflux = 10**(-0.4*(rmag-22.5))
        zflux = 10**(-0.4*((rmag-rzmag)-22.5))
        flux = np.vstack([gflux,rflux,zflux])
        
        # Randomly generate parameters.
        disk_r50 = np.random.uniform(min_disk_r50,max_disk_r50,nfake)
        bulge_r50 = np.random.uniform(min_bulge_r50,max_bulge_r50,nfake)
        nbulge = np.random.uniform(nbulgemin,nbulgemax,nfake)
        ndisk = np.random.uniform(ndiskmin,ndiskmax,nfake)
        bulge_frac = np.random.uniform(min_bulge_frac,max_bulge_frac,nfake)
        phi = np.random.uniform(phimin,phimax,nfake)
        axisratio = np.random.uniform(axismin,axismax,nfake)

        # Create a fits table containing the arrays of the randomly generated parameters.
        tbhdu = fits.BinTableHDU.from_columns([
            fits.Column(name='ra',format='f4',array=ra),
            fits.Column(name='dec',format='f4',array=dec),
            fits.Column(name='rmag',format='f4',array=rmag),
            fits.Column(name='grmag',format='f4',array=grmag),
            fits.Column(name='rzmag',format='f4',array=rzmag),
            fits.Column(name='nbulge',format='f4',array=nbulge),
            fits.Column(name='ndisk',format='f4',array=ndisk),
            fits.Column(name='disk_r50',format='f4',array=disk_r50),
            fits.Column(name='bulge_r50',format='f4',array=bulge_r50),
            fits.Column(name='bulge_frac',format='f4',array=bulge_frac),
            fits.Column(name='phi',format='f4',array=phi),
            fits.Column(name='axisratio',format='f4',array=axisratio)])
        tbhdu.writeto(fake_decals_dir+'decals_fake_priors.fits',clobber=True)

        # Get all the CCDs that touch this brick
        targetwcs = wcs_for_brick(decals.get_brick_by_name(thisbrick))
        allccds = decals.ccds_touching_wcs(targetwcs)
        # Testing (to be removed)
        allccds = allccds[:3]
        

        # Put all ccds in the proper directory; if the directory does not exist, create it.
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
                

        for ccd in allccds:
            print(ccd.ra,ccd.dec,ccd.filter,ccd.calname,ccd.cpimage)
            calname = ccd.calname.strip()
            filt = ccd.filter.strip()

            # Name the imags and place them in the proper directory.
            psffile = os.path.join(decals_dir,'calib','decam','psfex',calname+'.fits')
            wcsfile = os.path.join(decals_dir,'calib','decam','astrom-pv',calname+'.wcs.fits')

            #psf = psfex.PsfEx.fromFits(psffile)
            initpsf = psfex.PsfEx(psffile,imw,imh,ny=13,nx=7,
                                  psfClass=GaussianMixtureEllipsePSF,K=2)

            # Read the pre-existing decals image.
            imfile = os.path.join(fake_decals_dir,'images',ccd.cpimage.strip())
            print('Reading {}'.format(imfile))
            im = galsim.fits.read(imfile,hdu=ccd.ccdnum)
            hdr = galsim.fits.FitsHeader(imfile,hdu=ccd.ccdnum)
            #image, hdulist, fin = galsim.fits.readFile(imfile,hdu=ccd.ccdnum)
            #image = fits.getdata(imfile,hdu=ccd.ccdnum)
            #hdr = fits.getheader(imfile,hdu=ccd.ccdnum)
            #im = galsim.Image(image.data)
            wcs, origin = galsim.wcs.readFromFitsHeader(fits.getheader(wcsfile))

            #inverse variance array
            ivarfile = imfile.replace('ooi','oow')
            print('Reading inverse variance array{}'.format(ivarfile))
            invvar = galsim.fits.read(ivarfile,hdu=ccd.ccdnum)

            # Loop,which converts ra/dec to wcs
            for iobj in range(nfake):
                print(iobj)
                pos = wcs.posToImage(galsim.CelestialCoord(
                    ra[iobj]*galsim.degrees,dec[iobj]*galsim.degrees))
                pos = galsim.PositionD(1050.1214,2054.21341)
                xpos = int(pos.x)
                ypos = int(pos.y)
                dx = pos.x-xpos
                dy = pos.y-ypos

                localwcs = wcs.local(image_pos=pos)
                #localscale = np.sqrt(local.dudx*local.dvdy)

                # Ensure that the galaxy will be within the parameters of the ccd

                    #psf = galsim.Gaussian(flux=1.0, sigma=1.0)
                psfim = PsfEx.instantiateAt(initpsf,xpos,ypos)[5:-5,5:-5]
                psf = galsim.Image(psfim)
                psf = galsim.InterpolatedImage(psf,scale=1.0,flux=1.0)
                psf_centroid = psf.centroid()
                psf = psf.shift(-psf_centroid.x,-psf_centroid.y)
                    
                    # Creates the galaxies.

                    #bulge = galsim.Sersic(n=nbulge[iobj],half_light_radius=ndisk[iobj],gsparams=gsparams,flux=flux[iband,iobj])
                    #disk = galsim.Sersic(ndisk[iobj],scale_radius=disk_r50[iobj])
                    #stamp = bulge_frac[iobj] * bulge + (1-bulge_frac[iobj]) * disk
                gal = galsim.Sersic(ndisk[iobj],scale_radius=disk_r50[iobj],flux=1E8)
                gal = gal.shear(q=axisratio[iobj],beta=phi[iobj]*galsim.degrees)
                gal = gal.shift(dx=-100.0,dy=75.0)
                    #gal = gal.shift(dx=dx,dy=dy)
                gal = galsim.Convolve([gal,psf])

                stamp = gal.drawImage(wcs=localwcs)
                stamp.setCenter(xpos,ypos)

                    # Sets the bounds of the image.   
                bounds = stamp.bounds & im.bounds
                # if bounds is not within the image, crop stamp so it will fit in im and update bounds else
                if not xpos>0 and xpos<2046 and ypos>0 and ypos<4094:
                    #How do you slice the image at the edge?
                im[bounds] += stamp[bounds]

            # Writes the images to the output directory.
            outfile = os.path.join(fake_decals_dir,'images',ccd.cpimage.strip())
            print('Updating extension {} of image {}'.format(ccd.ccdnum,outfile))
            fits.update(imfile,im.array,ext=ccd.ccdnum,header=fits.Header(hdr.items()))

            galsim.fits.write(im,file_name='/home/desi3/junk.fits',clobber=True)

if __name__ == "__main__":
    decals_makefake()

# My to do list:
    # Work on pixel scale (balrog)
    # Offset?
    # PSF
   
    
