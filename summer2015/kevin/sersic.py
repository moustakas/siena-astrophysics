#!/usr/bin/env python

import galsim

def galaxy2():
    
    pixel_scale = 0.25
    half = 20.0
    n = 2.0
    flux = 1e5
    q = 0.2
    beta = 75 * galsim.degrees

    gal = galsim.Sersic(n=n, half_light_radius=half, flux=flux)
    psf = galsim.Gaussian(flux=1.0, sigma=1.0)
    noise = 0.001
    gal_shape = galsim.Shear(q=q, beta=beta)
    gal = gal.shear(gal_shape)
    image = galsim.Convolve([gal,psf])
    image = gal.drawImage(scale=pixel_scale)
    
    #image.addNoise(galsim.GaussianNoise(sigma=noise))
    print('Writing Test Galaxy')
    image.write('sersic1.fits',clobber=True)

if __name__ == "__main__":
    galaxy2()

