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
    nfake = 5
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

    stampwidth = 55 # postage stamp width [pixels, roughly 14 arcsec]
    
    # Set directories.
    decals_dir = os.getenv('DECALS_DIR')
    fake_decals_dir = os.getenv('FAKE_DECALS_DIR')

    brickinfo = fits.getdata(os.path.join(decals_dir,'decals-bricks.fits'),1)
    ccdinfo = fits.getdata(os.path.join(decals_dir,'decals-ccds.fits'),1)
    zpts =  fits.getdata(os.path.join(decals_dir,'calib','decam','photom','zeropoints.fits'))
    zptfilename = np.array(list(os.path.basename(ff).strip() for ff in zpts['FILENAME']))
    
    # Create an array for both the bricks and bands.
    mybricks = ['2428p117']
    #mybricks = ['2426p197','2428p117']
    band = np.array(['g','r','z'])

    # initialize some classes we'll need
    decals = Decals()
    
    #  Loop randomly generating priors.
    for thisbrick in mybricks:
        print('Getting info on brick {}'.format(thisbrick))
        info = brickinfo[np.where((brickinfo['brickname']==thisbrick)*1==1)]
        
        # Randomly generate location and flux parameters.
        ra = np.random.uniform(info['ra1'],info['ra2'],nfake)
        dec = np.random.uniform(info['dec1'],info['dec2'],nfake)

        print('Hack!')
        ra = np.random.uniform(242.68343197690362-1000*0.262/3600,242.68343197690362+1000*0.262/3600,nfake)
        dec = np.random.uniform(11.933249806710897-2000*0.262/3600,11.933249806710897+2000*0.262/3600,nfake)

        rmag = np.random.uniform(rmagmin,rmagmax,nfake)
        grmag = np.random.uniform(grmin,grmax,nfake)
        rzmag = np.random.uniform(rzmin,rzmax,nfake)

        # Calculate the g, r, and z band fluxes and stack them in an array.
        gflux = 10**(-0.4*(rmag+grmag))
        rflux = 10**(-0.4*rmag)
        zflux = 10**(-0.4*(rmag-rzmag))
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
        fakeinfo = fits.BinTableHDU.from_columns([
            fits.Column(name='ra',format='f4',array=ra),
            fits.Column(name='dec',format='f4',array=dec),
            fits.Column(name='rmag',format='f4',array=rmag),
            fits.Column(name='grmag',format='f4',array=grmag),
            fits.Column(name='rzmag',format='f4',array=rzmag),
            #fits.Column(name='flux',format='f4',array=flux),
            fits.Column(name='nbulge',format='f4',array=nbulge),
            fits.Column(name='ndisk',format='f4',array=ndisk),
            fits.Column(name='disk_r50',format='f4',array=disk_r50),
            fits.Column(name='bulge_r50',format='f4',array=bulge_r50),
            fits.Column(name='bulge_frac',format='f4',array=bulge_frac),
            fits.Column(name='phi',format='f4',array=phi),
            fits.Column(name='axisratio',format='f4',array=axisratio)])
        fakeinfo.writeto(fake_decals_dir+'decals_fake_priors.fits',clobber=True)

        # Get all the CCDs that touch this brick
        targetwcs = wcs_for_brick(decals.get_brick_by_name(thisbrick))
        allccds = decals.ccds_touching_wcs(targetwcs)
        # Testing (to be removed)
        allccds = allccds[:1]
        
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

            # Get filenames.
            cpname = os.path.basename(ccd.cpimage).strip().replace('.fz','')
            calname = ccd.calname.strip()
            thisband = np.where((ccd.filter==band)*1)

            imfile = os.path.join(fake_decals_dir,'images',ccd.cpimage.strip())
            ivarfile = imfile.replace('ooi','oow')
            psffile = os.path.join(decals_dir,'calib','decam','psfex',calname+'.fits')
            wcsfile = os.path.join(decals_dir,'calib','decam','astrom-pv',calname+'.wcs.fits')

            # Read the data.
            print('Reading {}'.format(imfile))
            im = galsim.fits.read(imfile,hdu=ccd.ccdnum)
            invvar = galsim.fits.read(ivarfile,hdu=ccd.ccdnum)
            hdr = galsim.fits.FitsHeader(imfile,hdu=ccd.ccdnum)

            # Get the WCS info and initialize the PSF
            wcs, origin = galsim.wcs.readFromFitsHeader(galsim.fits.FitsHeader(wcsfile))
            initpsf = psfex.PsfEx(psffile,im.xmax,im.ymax,ny=13,nx=7,
                                  psfClass=GaussianMixtureEllipsePSF,K=2)

            # Get the zeropoint info for this image and CCD
            zptinfo = zpts[np.where(((zptfilename==cpname)*1)*((ccd.extname==zpts['CCDNAME'])*1))]
            magzpt = zptinfo['ZPT'] + 2.5*np.log10(zptinfo['EXPTIME'])

            # Loop on each fake galaxy.
            for iobj in range(nfake):
                print(iobj, ra[iobj], dec[iobj])

                #ra[iobj] = ccd.ra-0.01
                #dec[iobj] = ccd.dec+0.01

                # Get the position of the galaxy on the CCD and build the PSF.
                pos = wcs.toImage(galsim.CelestialCoord(
                    ra[iobj]*galsim.degrees,dec[iobj]*galsim.degrees))
                xpos = int(pos.x)
                ypos = int(pos.y)

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
                gal = galsim.Sersic(ndisk[iobj],scale_radius=disk_r50[iobj],
                                    flux=1000*float(flux[thisband,iobj]*10**(-0.4*magzpt)))
                gal = gal.shear(q=axisratio[iobj],beta=phi[iobj]*galsim.degrees)
                gal = galsim.Convolve([gal,psf])
                #gal = gal.shift(pos.x-xpos,pos.y-ypos) # apply sub-pixel shift
                
                stamp = gal.drawImage(wcs=wcslocal,method='no_pixel')
                #stamp = gal.drawImage(nx=stampwidth,ny=stampwidth,wcs=wcslocal,
                #                      offset=offset,method='no_pixel')
                stamp.setCenter(xpos,ypos)

                overlap = stamp.bounds & im.bounds
                if (overlap.xmax>=0 and overlap.ymax>=0 and overlap.xmin<=im.xmax and
                    overlap.ymin<=im.ymax and overlap.area()>0):
                    print('Adding object {} to ...'.format(iobj))
                    print(xpos, ypos)
                    im[overlap] += stamp[overlap]

            # Writes the images to the output directory.
            outfile = os.path.join(fake_decals_dir,'images',ccd.cpimage.strip())
            print('Updating extension {} of image {}'.format(ccd.ccdnum,outfile))
            fits.update(imfile,im.array,ext=ccd.ccdnum,header=fits.Header(hdr.items()))

            galsim.fits.write(im,file_name='junk.fits',clobber=True)

if __name__ == "__main__":
    decals_makefake()

# My to do list:
    # Work on pixel scale (balrog)
    # Offset?
    # PSF
   
    
