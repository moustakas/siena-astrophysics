#!/usr/bin.env python

import galsim

def DeVaucouleur():
    r50 = 4.0
    flux = 100.0
    trunc = 10
    scale = 0.25
    noise = 10
    #psf = galsim.Gaussian(flux=1.0,sigma=1.0)
    gal = galsim.DeVaucouleurs(half_light_radius=r50,flux=flux,trunc=trunc)
    image = gal.drawImage(scale=scale)
    #image = galsim.Convolve([gal,psf])
    #image.addNoise(galsim.GaussianNoise(sigma=noise))
    image.write('devau1.fits')

if __name__ == "__main__":
    DeVaucouleur()
