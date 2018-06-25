#!/usr/bin/env python

import galsim

def blank():
    
    pixel_scale = 0.262
    gal = galsim.Gaussian(half_light_radius=10 ,flux=0)
    noise = 0.01
    psf = 0.001
    image = gal.drawImage(scale=pixel_scale)    
    image.addNoise(galsim.GaussianNoise(sigma=noise))
 
    print('Writing Test Galaxy')
    image.write('blank.fits',clobber=True)

if __name__ == "__main__":
    blank()
