#!/usr/bin/env python

import galsim
import os
def bulge():
    
    gal_flux = 1.e6        
    bulge_n = 3.5       
    bulge_re = 2.3      
    disk_n = 1.5        
    disk_r0 = 0.85      
    bulge_frac = 0.3    
    gal_q = 0.73        
    gal_beta = 23        
    atmos_fwhm=2.1       
    atmos_e = 0.13       
    atmos_beta = 0.81    
    opt_defocus=0.53     
    opt_a1=-0.29         
    opt_a2=0.12          
    opt_c1=0.64          
    opt_c2=-0.33         
    opt_obscuration=0.3  
    lam = 800            
    tel_diam = 4.        
    pixel_scale = 0.23   
    image_size = 64      
    wcs_g1 = -0.02       
    wcs_g2 = 0.01        
    sky_level = 2.5e4    
    gain = 1.7           
    read_noise = 0.3
    random_seed = 1314662  

    rng = galsim.BaseDeviate(random_seed) #random number generator
    bulge = galsim.Sersic(bulge_n, half_light_radius=bulge_re)
    disk = galsim.Sersic(disk_n, scale_radius=disk_r0)
    gal = bulge_frac * bulge + (1-bulge_frac) * disk
    gal = gal.withFlux(gal_flux)
    gal_shape = galsim.Shear(q=gal_q, beta=gal_beta*galsim.degrees)
    gal = gal.shear(gal_shape)   
    atmos = galsim.Kolmogorov(fwhm=atmos_fwhm)
    atmos = atmos.shear(e=atmos_e, beta=atmos_beta*galsim.radians)
    lam_over_diam = lam * 1.e-9 / tel_diam # radians
    lam_over_diam *= 206265  # arcsec
    optics = galsim.OpticalPSF(lam_over_diam, 
                               defocus = opt_defocus,
                               coma1 = opt_c1, coma2 = opt_c2,
                               astig1 = opt_a1, astig2 = opt_a2,
                               obscuration = opt_obscuration) 
    wcs = galsim.ShearWCS(scale=pixel_scale, shear=galsim.Shear(g1=wcs_g1, g2=wcs_g2))   
    psf = galsim.Convolve([atmos, optics])
    final = galsim.Convolve([psf, gal])
    image = galsim.ImageF(image_size, image_size)
    image_epsf = galsim.ImageF(image_size, image_size)
    psf.drawImage(image_epsf, wcs=wcs)
    image_opticalpsf = optics.drawImage(method='sb')
    image += sky_level * pixel_scale**2
    image.addNoise(galsim.CCDNoise(rng, gain=gain, read_noise=read_noise))
    image -= sky_level * pixel_scale**2
    file_name = os.path.join('output', 'demo3.fits')
    file_name_epsf = os.path.join('output','demo3_epsf.fits')
    file_name_opticalpsf = os.path.join('output','demo3_opticalpsf.fits')
    image.write(file_name)
    image_epsf.write(file_name_epsf)
    image_opticalpsf.write(file_name_opticalpsf)

if __name__ == "__main__":
    bulge()
