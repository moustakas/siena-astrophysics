#!/usr/bin/env python

import galsim
import numpy as np
import os
import math
def InsertSimulatedGalaxies():

    nx_tiles = 5                   
    ny_tiles = 5                   
    stamp_xsize = 40               
    stamp_ysize = 40               

    random_seed = 6424512          

    pixel_scale = 1.0               # arcsec / pixel
    sky_level = 1.e6                # ADU / arcsec^2

    if not os.path.isdir('output'):
        os.mkdir('output')
    psf_file_name = os.path.join('output','g08_psf.fits')
    psf_beta = 3                 
    psf_fwhm = 2.85                 # arcsec (=pixels)
    psf_trunc = 2.*psf_fwhm         # arcsec (=pixels)
    psf_e1 = -0.019              
    psf_e2 = -0.007              

    gal_file_name = os.path.join('output','g08_gal.fits')
    gal_signal_to_noise = 200      
    gal_resolution = 0.98           # r_gal / r_psf (use r = half_light_radius)
    gal_ellip_rms = 0.2             # using "distortion" definition of ellipticity:
                                    #   e = (a^2-b^2)/(a^2+b^2), where a and b are the 
                                    #   semi-major and semi-minor axes, respectively.
    gal_ellip_max = 0.6             # Maximum value of e, to avoid getting near e=1.
    gal_g1 = 0.5                  # Applied shear, using normal shear definition:
    gal_g2 = -0.008                 #   g = (a-b)/(a+b)

    shift_radius = 1.0              # arcsec (=pixels)

  
    # Define the PSF profile
    psf = galsim.Moffat(beta=psf_beta, fwhm=psf_fwhm, trunc=psf_trunc)

    # When something can be constructed from multiple sizes, e.g. Moffat, then
    # you can get any size out even if it wasn't the way the object was constructed.
    # In this case, we extract the half-light radius, even though we built it with fwhm.
    # We'll use this later to set the galaxy's half-light radius in terms of a resolution.
    psf_re = psf.getHalfLightRadius()

    psf = psf.shear(e1=psf_e1,e2=psf_e2)
  

    # Define the galaxy profile

    # First figure out the size we need from the resolution
    gal_re = psf_re * gal_resolution

    # Make the galaxy profile starting with flux = 1.
    gal = galsim.Exponential(flux=1., half_light_radius=gal_re)
  

    # This profile is placed with different orientations and noise realizations
    # at each postage stamp in the gal image.
    gal_image = galsim.ImageF(stamp_xsize * nx_tiles-1 , stamp_ysize * ny_tiles-1,
                              scale=pixel_scale)
    psf_image = galsim.ImageF(stamp_xsize * nx_tiles-1 , stamp_ysize * ny_tiles-1,
                              scale=pixel_scale)

    shift_radius_sq = shift_radius**2

    first_in_pair = True  # Make pairs that are rotated by 90 degrees

    k = 0
    for iy in range(ny_tiles):
        for ix in range(nx_tiles):
            ud = galsim.UniformDeviate(random_seed+k)

            gd = galsim.GaussianDeviate(ud, sigma=gal_ellip_rms)

            b = galsim.BoundsI(ix*stamp_xsize+1 , (ix+1)*stamp_xsize-1, 
                               iy*stamp_ysize+1 , (iy+1)*stamp_ysize-1)
            sub_gal_image = gal_image[b]
            sub_psf_image = psf_image[b]

            # Great08 randomized the locations of the two galaxies in each pair,
            # but for simplicity, we just do them in sequential postage stamps.
            if first_in_pair:
                # Use a random orientation:
                beta = ud() * 2. * math.pi * galsim.radians

                # Determine the ellipticity to use for this galaxy.
                ellip = 1
                while (ellip > gal_ellip_max):
                    # Don't do `ellip = math.fabs(gd())`
                    # Python basically implements this as a macro, so gd() is called twice!
                    val = gd()
                    ellip = math.fabs(val)

                first_in_pair = False
            else:
                # Use the previous ellip and beta + 90 degrees
                beta += 90 * galsim.degrees
                first_in_pair = True

            # Make a new copy of the galaxy with an applied e1/e2-type distortion 
            # by specifying the ellipticity and a real-space position angle
            this_gal = gal.shear(e=ellip, beta=beta)

            # Apply the gravitational reduced shear by specifying g1/g2
            this_gal = this_gal.shear(g1=gal_g1, g2=gal_g2)

            # Apply a random shift_radius:
            rsq = 2 * shift_radius_sq
            while (rsq > shift_radius_sq):
                dx = (2*ud()-1) * shift_radius
                dy = (2*ud()-1) * shift_radius
                rsq = dx**2 + dy**2

            this_gal = this_gal.shift(dx,dy)
            # Note that the shifted psf that we create here is purely for the purpose of being able
            # to draw a separate, shifted psf image.  We do not use it when convolving the galaxy
            # with the psf.
            this_psf = psf.shift(dx,dy)

            # Make the final image, convolving with the (unshifted) psf
            final_gal = galsim.Convolve([psf,this_gal])

            # Draw the image
            final_gal.drawImage(sub_gal_image)

            sky_level_pixel = sky_level * pixel_scale**2
            noise = galsim.PoissonNoise(ud, sky_level=sky_level_pixel)
            sub_gal_image.addNoiseSNR(noise, gal_signal_to_noise)

            # Draw the PSF image
            # No noise on PSF images.  Just draw it as is.
            this_psf.drawImage(sub_psf_image)

            # For first instance, measure moments
            if ix==0 and iy==0:
                psf_shape = sub_psf_image.FindAdaptiveMom()
                temp_e = psf_shape.observed_shape.e
                if temp_e > 0.0:
                    g_to_e = psf_shape.observed_shape.g / temp_e
                else:
                    g_to_e = 0.0
         
            x = b.center().x
            y = b.center().y
            k = k+1

  

    # Now write the images to disk.
    psf_image.write(psf_file_name)
  

    gal_image.write(gal_file_name)



if __name__ == "__main__":
    InsertSimulatedGalaxies()
    
