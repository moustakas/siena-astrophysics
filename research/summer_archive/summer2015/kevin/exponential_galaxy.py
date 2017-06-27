#!/usr/bin/env python

#Creating a postage stamp of an exponential galaxy
#If we're lucky, import it into a coordinate in an image

#Import stuff
import sys
import os

#Define program
def galaxy(gauss=False,exp=True):
    import galsim
    psf_sigma = 1.
    psf = galsim.Gaussian(flux=1., sigma=psf_sigma)
    pixel_scale = 0.25
    noise = 0.01   
    
    if gauss:
        #Insert galaxy parameters
        gal_flux = 1.e5
        gal_sigma = 2.
        gal = galsim.Gaussian(flux=gal_flux, sigma=gal_sigma)

    if exp:
        half_radius = 10.0 #Arcseconds
        flux = 1000.0
        gal = galsim.Exponential(half_light_radius=half_radius,flux=flux)

    image = galsim.Convolve([gal,psf])
    image = gal.drawImage(scale=pixel_scale)
    #image.addNoise(galsim.GaussianNoise(sigma=noise))
    print('Writing Test Galaxy')
    image.write('test.fits',clobber=True)

if __name__ == "__main__":
    galaxy()
