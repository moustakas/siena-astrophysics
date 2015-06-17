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

    gsparams = galsim.GSParams(maximum_fft_size=2L**30L)

    # Set parameter ranges.
    nfake = 20
    min_disk_r50 = 0.1
    max_disk_r50 = 2.0
    min_bulge_r50 = 0.1
    max_bulge_r50 = 1.0
    nbulgemin = 3.0
    nbulgemax = 5.0
    ndiskmin = 1.0
    ndiskmax = 1.0
    min_bulge_frac = 0.0
    max_bulge_frac = 1.0
    phimin = 0.0
    phimax = 180.0
    axismin = 0.2
    axismax = 1.0
    rmagmin = 17.0
    rmagmax = 19.0
    grmin = -0.3
    grmax = 0.5
    rzmin = 0.0
    rzmax = 1.5

    avgpixscale = 0.262 # [arcsec/pixel]
    stampwidth = 95 # postage stamp width [pixels, roughly 14 arcsec]
    stampbounds = galsim.BoundsD(0,stampwidth,0,stampwidth)
    
    # Set directories.
    decals_dir = os.getenv('DECALS_DIR')
    fake_decals_dir = os.getenv('FAKE_DECALS_DIR')

    brickinfo = fits.getdata(os.path.join(decals_dir,'decals-bricks.fits'),1)
    ccdinfo = fits.getdata(os.path.join(decals_dir,'decals-ccds.fits'),1)
    zpts =  fits.getdata(os.path.join(decals_dir,'decals-ccds-zeropoints.fits'))
    
    # Create an array for both the bricks and bands.
    mybricks = ['2428p117']
    #mybricks = ['2426p197','2428p117']
    band = np.array(['g','r','z'])

    # initialize some classes we'll need
    decals = Decals()
    
    #  Loop randomly generating priors.
    for thisbrick in mybricks:
        print('Working on brick {}'.format(thisbrick))
        info = brickinfo[np.where((brickinfo['brickname']==thisbrick)*1==1)]
        
        # Randomly generate location and flux parameters.
        ramin = info['ra1']
        ramax = info['ra2']
        decmin = info['dec1']
        decmax = info['dec2']

#       ramin = 242.68343197690362-0.1
#       ramax = 242.68343197690362+0.1
#       decmin = 11.933249806710897-0.1
#       decmax = 11.933249806710897+0.1

        ## Temporarily focus on the center of the brick.
        ramin = info['ra']-500*avgpixscale/3600
        ramax = info['ra']+500*avgpixscale/3600
        decmin = info['dec']-500*avgpixscale/3600
        decmax = info['dec']+500*avgpixscale/3600

        ra = np.random.uniform(ramin,ramax,nfake)
        dec = np.random.uniform(decmin,decmax,nfake)

        rmag = np.random.uniform(rmagmin,rmagmax,nfake)
        grmag = np.random.uniform(grmin,grmax,nfake)
        rzmag = np.random.uniform(rzmin,rzmax,nfake)

        # Calculate the g, r, and z band fluxes and stack them in an array.
        gflux = 10**(-0.4*(rmag+grmag)) # [nanomaggies]
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
        #allccds = allccds[:1]
        
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
            for iobj in range(nfake):
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

if __name__ == "__main__":
    decals_makefake()

# My to do list:
    # Work on pixel scale (balrog)
    # Offset?
    # PSF
   
    
