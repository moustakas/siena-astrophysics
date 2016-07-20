#!/usr/bin/env python

"""Generate the CUTE input files from the SDSS catalogs and randoms posted at
http://data.sdss3.org/sas/dr11/boss/lss/

"""
from __future__ import division, print_function

import os
import sys
import argparse
import glob
import logging as log

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import WMAP7
from matplotlib.colors import LogNorm

def plotmqh(monopole,quadrupole,hexadecapole,rrange):
    plt.plot(rrange, monopole*rrange**2, 'ko')

def compute_monopole(mu, r, xirm):
    xirm = xirm*1.0
    Bxirm = np.reshape(xirm,[40,50]) # generalize
    xr = 0.025 # find out about this factor (factor of 1/40)
    monopole = xr*np.trapz(Bxirm, axis=0) # go through the math of why
    return monopole

def compute_quadrupole(mu, r, xirm):
    xirm = xirm*(3*(mu*mu)-1.0)*(5/2)
    Bxirm = np.reshape(xirm,[40,50])
    xr = 0.025
    quadrupole = xr*np.trapz(Bxirm)
    return quadrupole

def compute_hexadecapole(mu, r, xirm):
    xirm = xirm*(35*(mu*mu*mu*mu)-(30*mu*mu)+3)/8
    Bxirm = np.reshape(xirm,[40,50])
    xr = 0.025
    hexadecapole = xr*np.trapz(Bxirm)
    return hexadecapole

def calc_fkp_weights(z, zmin, zmax): 
    NRB = 200 
    NGC = 6308
    SGC = 2069
    SURVEY_SIZE = NGC
    FULL_AREA = 41253
    dz = zmax - zmin 
    red_interval = dz/NRB
    red_markers = []
    red_vol = []
    bin_num = []
    bin_sum = []
    wfkp = []

    for ii in range(NRB+1):
        red_markers.append(zmin + ii*red_interval) 
        if ii >= 1:
            red_vol.append((4/3)*np.pi*((WMAP7.comoving_distance(red_markers[ii]).value*0.704)**3-
                                        (WMAP7.comoving_distance(red_markers[ii-1]).value*0.704)**3)
                                        *(SURVEY_SIZE/FULL_AREA))
    for ii in range(len(z)):
        bin_num.append(int(np.floor(NRB * (z[ii]-zmin)/dz))) 
    bin_num = np.asarray(bin_num)

    for ii in range(NRB):
        bin_sum.append(len(np.where(bin_num==ii)[0])) 

    bin_sum = np.asarray(bin_sum)
    red_vol = np.asarray(red_vol)
    nbar_slice = bin_sum/red_vol

    for ii in range(len(z)):
        wfkp.append(1/(1+20000*nbar_slice[bin_num[ii]]))
                 
    return wfkp

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--dr', type=str, default='dr11', help='Specify the SDSS data release.')
    parser.add_argument('--parse', action='store_true', help='Parse the input datafiles (do just once).')
    parser.add_argument('--docute', action='store_true', help='Run CUTE.')
    parser.add_argument('--qaplots', action='store_true', help='Generate QAplots.')
    parser.add_argument('--corrtype', type=str, default='monopole', help='Specify correlation type (monopole|3D_ps|3D_rm).')
    parser.add_argument('--cosmo', type=int, default=1, help='Adopted cosmology (1|2)')

    args = parser.parse_args()

    key = 'LSS_BOSS'
    if key not in os.environ:
        log.fatal('Required ${} environment variable not set'.format(key))
        return 0

    CUTEdir = os.path.join(os.getenv('CUTE'))
    drdir = os.path.join(os.getenv('LSS_BOSS'), args.dr)
    randomsdir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'randoms')
    datafile = os.path.join(drdir, args.dr+'_cmass.dat')
    randomfile = os.path.join(drdir, 'parsed', '{}_cmass_random'.format(args.dr))
    baseoutfile = os.path.join(drdir, 'cuteout', args.corrtype, 'dr11_2pt_{}'.format(args.corrtype))
    paramfile = os.path.join(drdir, 'param', 'dr11_{}_'.format(args.corrtype))
    randomslist = glob.glob(os.path.join(randomsdir, '*.dat'))

    if args.corrtype == '3D_ps':
        npibins = 150
        nsigbins = 150
        maxpi = 150
        maxsig = 150

    # Parse the input data and write out CUTE-compatible files.
    if args.parse:
        allspecz = fits.getdata(os.path.join(drdir, 'galaxy_DR11v1_CMASS_North.fits.gz'), 1)
        keep = np.where((allspecz['Z']>0.43)*(allspecz['Z']<0.7))[0]
        specz = allspecz[keep]
        ngal = len(keep)
        datafkp = calc_fkp_weights(specz['Z'], 0.43, 0.7)
        data = np.zeros((ngal, 4))
        data[:,0] = specz['RA']
        data[:,1] = specz['DEC']
        data[:,2] = specz['Z']
        data[:,3] = datafkp*specz['WEIGHT_SYSTOT']*(specz['WEIGHT_NOZ']+specz['WEIGHT_CP']-1)
        print('Writing {}'.format(datafile))
        log.info('Writing {}'.format(datafile))
        np.savetxt(datafile, data)
	
        for item in range(len(randomslist)):
            ra, dec, z, ipoly, wboss, wcp, wzf, veto = \
              np.loadtxt(os.path.join(randomsdir, randomslist[item]), unpack=True)
            keep = np.where(veto==1)[0]
            nobj = len(keep)
            rand = np.zeros((nobj,4))
            rand[:,0] = ra[keep]
            rand[:,1] = dec[keep]
            rand[:,2] = z[keep]
            randfkp = calc_fkp_weights(rand[:,2], 0.43, 0.7)
            rand[:,3] = randfkp*(wcp[keep]+wzf[keep]-1)
            log.info('Writing {}'.format(randomfile+'_'+args.corrtype+'_fkp_{}.dat'.format(item+4001)))
            print('Writing {}'.format(randomfile+'_'+args.corrtype+'_fkp_{}.dat'.format(item+4001)))
            np.savetxt(randomfile+'_fkp_{}.dat'.format(item+4001), rand) # rename all of the randomefile+3D_rm files
                      
    if args.docute:

        # Select the cosmological parameters
        if args.cosmo == 1:
            omega_M = 0.3
            omega_L = 0.7
            ww = -1
            
        if args.cosmo == 2:
            omega_M = 0.274
            omega_L = 1.0 - omega_M
            ww = -1
    
        for item in range(len(randomslist)):
            # Create a unique filename for each parameeter file
            newfile = paramfile+'fkp_{}.param'.format(item+4001)

            # Write the parameter file; constants, and then conditionals
            pfile = open(newfile, 'w')
            pfile.write('data_filename= '+datafile+'\n')
            pfile.write('random_filename= '+randomfile+'_fkp_{}.dat'.format(item+4001)+'\n') # get rid of 3D_rm name
            pfile.write('mask_filename= junk\n')
            pfile.write('z_dist_filename= junk\n')
            pfile.write('output_filename= '+baseoutfile+'_fkp_{}.dat'.format(item+4001)+'\n')
            pfile.write('corr_type= '+args.corrtype+'\n')
            pfile.write('num_lines= all\n')
            pfile.write('corr_estimator= LS\n')
            
            if args.corrtype == 'monopole':
                pfile.write('input_format= 2\n')
                pfile.write('np_rand_fact= 8\n')
                pfile.write('omega_M= {}\n'.format(omega_M))
                pfile.write('omega_L= {}\n'.format(omega_L))
                pfile.write('w= {}\n'.format(ww))
                pfile.write('log_bin= 1\n') # changed from 0 to 1
                pfile.write('n_logint= 23\n')
                pfile.write('dim1_max= 150\n')
                pfile.write('dim1_nbin= 75\n')
                pfile.write('dim2_max= 150\n')
                pfile.write('dim2_nbin= 75\n')
                pfile.write('dim3_min= 0.4\n')
                pfile.write('dim3_max= 0.7\n')
                pfile.write('dim3_nbin= 1\n')
                pfile.write('radial_aperture= 1\n')
                pfile.write('use_pm= 1\n')
                pfile.write('n_pix_sph= 2048\n')
                
            if args.corrtype == '3D_rm':
                pfile.write('input_format= 2\n')
                pfile.write('np_rand_fact= 9.5217\n')
                pfile.write('omega_M= {}\n'.format(omega_M))
                pfile.write('omega_L= {}\n'.format(omega_L))
                pfile.write('w= {}\n'.format(ww))
                pfile.write('log_bin= 1\n') # changed from 0 to 1
                pfile.write('n_logint= 0\n')
                pfile.write('dim1_max= 160\n')
                pfile.write('dim1_nbin= 50\n')
                pfile.write('dim2_max= 1\n')
                pfile.write('dim2_nbin= 200\n')
                pfile.write('dim3_min= 0.4\n')
                pfile.write('dim3_max= 0.7\n')
                pfile.write('dim3_nbin= 1\n')
                pfile.write('radial_aperture= 1\n')
                pfile.write('use_pm= 0\n')
                pfile.write('n_pix_sph= 2048\n')

            if args.corrtype == '3D_ps':
                pfile.write('input_format= 2\n')
                pfile.write('np_rand_fact= 9.5217\n')
                pfile.write('omega_M= {}\n'.format(omega_M))
                pfile.write('omega_L= {}\n'.format(omega_L))
                pfile.write('w= {}\n'.format(ww))
                pfile.write('log_bin= 0\n')
                pfile.write('n_logint= 0\n')
                pfile.write('dim1_max= {}\n'.format(maxsig)) # Maximum Radial Separation
                pfile.write('dim1_nbin= {}\n'.format(nsigbins))
                pfile.write('dim2_max= {}\n'.format(maxpi)) # Maximum Transverse Separation
                pfile.write('dim2_nbin= {}\n'.format(npibins))
                pfile.write('dim3_min= 0.4\n')
                pfile.write('dim3_max= 0.7\n')
                pfile.write('dim3_nbin= 1\n')
                pfile.write('radial_aperture= 1\n')
                pfile.write('use_pm= 0\n')
                pfile.write('n_pix_sph= 2048\n')

            pfile.close()
            # Run CUTE, passing the newly created parameter file
            os.system('CUTE '+newfile)

    if args.qaplots:
        
        anderson1 = os.path.join(drdir, 'Anderson_2013_CMASSDR11_corrfunction_x0x2_prerecon.dat')
        and_rad,and_mono,and_quad = np.loadtxt(anderson1, unpack=True)

        if args.corrtype == 'monopole':
            for item in range(len(randomslist)):
                thisout = outfile+'{}.dat'.format(item+4001)
                rad, xi, xierr, DD, DR, RR = np.loadtxt(thisout, unpack=True)
                plt.scatter(rad, xi*rad**2)
                #plt.axis([-5, 155, 0, 120])
                plt.xlabel('$\mathrm{\ r \ (Mpc)}$')
                plt.ylabel(r'$\mathrm{\ r^2 * \xi}$')
                #plt.savefig(os.path.join(drdir,'qaplots','power_spectrum_monopole_{}.png'.format(item+4001)))
                plt.show()
                
        if args.corrtype == '3D_rm':
            for item in range(len(randomslist)):
                thisout = outfile+'fkp_{}.dat'.format(item+4001)
                mu, rad, xi, xierr, DD, DR, RR = np.loadtxt(thisout, unpack=True)
                rad = np.linspace(2, 198, 40)
                # rad = np.linspace(2, 198, 40)
                # rad = rad.reshape((50,40))
                monopole = compute_monopole(mu, rad, xi)
                quadrupole = compute_quadrupole(mu, rad, xi)
                plotmqh(monopole,quadrupole,hex1,rad)
                plt.xlabel('$\mathrm{\ r \ (Mpc \,  h^{-1})}$')
                plt.ylabel(r'$\mathrm{\ r^2 \xi(r)}$')
            plt.plot(and_rad, (and_mono)*and_rad**2, 'r-')
            plt.show()
                
        if args.corrtype == '3D_ps':
            xi = np.zeros((len(randomslist), nsigbins, npibins))
            for item in range(len(randomslist)):
                thisout = outfile+'fkp_{}.dat'.format(item+4001)
                pi, sigma, thisxi, xierr, DD, DR, RR = np.loadtxt(thisout, unpack=True)
                xi = thisxi.reshape(nsigbins, npibins)
            #xi = np.mean(xi)
            bigxi = np.zeros((nsigbins*2, npibins*2))
            xi2d = xi.reshape(nsigbins, npibins)
            bigxi[:nsigbins, :npibins] = np.fliplr(np.flipud(xi2d))
            bigxi[nsigbins:2*nsigbins, :npibins] = np.fliplr(xi2d)
            bigxi[:nsigbins, npibins:2*npibins] = np.flipud(xi2d)
            bigxi[nsigbins:2*nsigbins, npibins:2*npibins] = xi2d
            pi2d = np.tile(np.arange(-(2*npibins-1), 2*npibins, 2), (1, 2*npibins)).reshape(2*npibins, 2*npibins)
            sig2d = np.rot90(pi2d)
            plt.pcolor(sig2d, pi2d, bigxi, norm=LogNorm()) ; plt.colorbar() ; plt.show()      
                            
if __name__ == "__main__":
    main()
