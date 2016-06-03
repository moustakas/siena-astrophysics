#!/usr/bin/env python

"""Generate the CUTE input files from the SDSS catalogs and randoms posted at
http://data.sdss3.org/sas/dr11/boss/lss/

"""
# monopole vs published
# calculate random weights and reimpliment data weights
# figure out the image
# covariance matrix

from __future__ import division, print_function

import os
import sys
import argparse
import glob
import logging as log
import numpy as np
import PIL
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import Planck13

def plotmqh(mono1,q1,hx1,rrange):
    plt.figure()
    plt.subplot(311)
    plt.plot(rrange, mono1*rrange**2, 'ko')
    plt.subplot(312)
    plt.plot(rrange, q1*rrange**2, 'ko')
    plt.subplot(313)
    plt.plot(rrange, hx1*rrange**2, 'ko')
    plt.show()

def compute_monopole(mu, r, xirm):
    xirm = xirm*1.0
    Bxirm = np.reshape(xirm,[40,50])
    xr = 0.025
    mono1 = xr*np.trapz(Bxirm)
    print(len(mono1))
    return mono1

def compute_quadrupole(mu, r, xirm):
    xirm = xirm*(3*(mu*mu)-1.0)*(5/2)
    Bxirm = np.reshape(xirm,[40,50])
    xr = 0.025
    q1 = xr*np.trapz(Bxirm)
    return q1

def compute_hexadecapole(mu, r, xirm):
    xirm = xirm*(35*(mu*mu*mu*mu)-(30*mu*mu)+3)/8
    Bxirm = np.reshape(xirm,[40,50])
    xr = 0.025
    hx1 = xr*np.trapz(Bxirm)
    return hx1

def calc_fkp_weights(z, zmax, zmin): # z is a list of redshifts, zmax and zmin are constants
    NRB = 200 # Number of redshift 
    SURVEY_SIZE = 8500 # dr11 survey size, square degrees
    FULL_AREA = 41253 # full area of the sky, square degrees
    dz = zmax - zmin # difference in redshift
    red_interval = dz/NRB # redshift width of each slice
    for ii in range(NRB+1):
        red_markers = zmin + ii*red_interval # slice
        if ii >= 1:
            red_vol = (Planck13.comoving_volume(red_markers[ii])-Planck13.comoving_volume(red_markers[ii-1]))*(SURVEY_SIZE/FULL_AREA) # find the volme of each slice
    for ii in z:
        bin_num = NRB * (z-zmin)/dz # assigns a bin number to each galaxy
        bin_sum = sum(np.where(bin_num==ii)) # finds the number of galaxies in each bin
    nbar_slice = bin_sum/red_vol # finds the mean number density of galaxies in a redslice
    wfkp = 1/(1+20000*nbar_slice) # finds the fkp weight of each source
    return weight

def covariance(rad, xi):
    cov = np.cov(rad)
    return blah

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--dr', type=str, default='dr11', help='Specify the SDSS data release.')
    parser.add_argument('--parse', action='store_true', help='Parse the input datafiles.')
    parser.add_argument('--docute', type=str, default='monopole', help='Run CUTE.')
    parser.add_argument('--qaplots', type=str, default=None, help='Generate QAplots.')
    parser.add_argument('--fkp', action='store_true', help='Calculate the fkp weights of the randoms.')

    args = parser.parse_args()

    key = 'LSS_BOSS'
    if key not in os.environ:
        log.fatal('Required ${} environment variable not set'.format(key))
        return 0

    CUTEdir = os.path.join(os.getenv('CUTE'))
    drdir = os.path.join(os.getenv('LSS_BOSS'), args.dr)
    randomsdir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'randoms')
    datafile = os.path.join(drdir, args.dr+'_cmass.dat')
    randomfile = os.path.join(drdir, 'parsed', args.dr+'_cmass_random')
    outfile = os.path.join(drdir, 'cuteout', 'dr11_2pt_'+args.docute+'_')
    paramfile = os.path.join(drdir, 'param', 'dr11_'+args.docute+'_')
    randomslist = glob.glob(os.path.join(randomsdir, '*.dat'))

    if args.fkp:
        calc_fkp_weights()
        
    # Parse the input data and write out CUTE-compatible files.
    if args.parse:

        allspecz = fits.getdata(os.path.join(drdir, 'galaxy_DR11v1_CMASS_North.fits.gz'), 1)
        keep = np.where((allspecz['Z']>0.43)*(allspecz['Z']<0.7))[0]
        specz = allspecz[keep]
        ngal = len(keep)

        data = np.zeros((ngal,4))
        data[:,0] = specz['RA']
        data[:,1] = specz['DEC']
        data[:,2] = specz['Z']
        data[:,3] = specz['WEIGHT_SYSTOT']*(specz['WEIGHT_NOZ']+specz['WEIGHT_CP']-1)
        # specz['WEIGHT_FKP']*specz['WEIGHT_SYSTOT']*(specz['WEIGHT_NOZ']+specz['WEIGHT_CP']-1)
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
            calc_fkp_weights(rand[:,2], 0.43, 0.7)
            rand[:,3] = wcp[keep]+wzf[keep]-1
            #log.info('Writing {}'.format(randomfile))
            print('Writing file {} of 4600'.format(item+4001))
            np.savetxt(randomfile+'{}.dat'.format(item+4001), rand)
          
    if args.docute == True and args.qaplots != None:
        for item in range(len(randomslist)):
            newfile = paramfile+'{}.param'.format(item+4001)

            pfile = open(newfile, 'w')
            pfile.write('data_filename= '+datafile+'\n')
            pfile.write('random_filename= '+randomfile+'{}.dat'.format(item+4001)+'\n')
            pfile.write('mask_filename= junk\n')
            pfile.write('z_dist_filename= junk\n')
            pfile.write('output_filename= '+outfile+'{}.dat'.format(item+4001)+'\n')
            pfile.write('corr_type= '+args.docute+'\n')
            pfile.write('num_lines= all\n')
            pfile.write('corr_estimator= LS\n')
            
            if args.docute == 'monopole':
                pfile.write('input_format= 2\n')
                pfile.write('np_rand_fact= 8\n')
                pfile.write('omega_M= 0.3\n')
                pfile.write('omega_L= 0.7\n')
                pfile.write('w= -1\n')
                pfile.write('log_bin= 0\n')
                pfile.write('n_logint= 10\n')
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
                
            if (args.docute == '3D_rm' or args.docute == '3D_ps'):
                pfile.write('input_format= 2\n')
                pfile.write('np_rand_fact= 9.5217\n')
                pfile.write('omega_M= 0.3\n')
                pfile.write('omega_L= 0.7\n')
                pfile.write('w= -1\n')
                pfile.write('log_bin= 0\n')
                pfile.write('n_logint= 0\n')
                pfile.write('dim1_max= 200\n')
                pfile.write('dim1_nbin= 50\n')
                pfile.write('dim2_max= 1\n')
                pfile.write('dim2_nbin= 40\n')
                pfile.write('dim3_min= 0.4\n')
                pfile.write('dim3_max= 0.7\n')
                pfile.write('dim3_nbin= 1\n')
                pfile.write('radial_aperture= 1\n')
                pfile.write('use_pm= 0\n')
                pfile.write('n_pix_sph= 2048\n')

            pfile.close()
            os.system('CUTE '+newfile)

    if args.qaplots:
            # Make rockin' plots and write out.
            if args.qaplots == 'monopole':
                for item in range(len(randomslist)):
                    thisout = outfile+'{}.dat'.format(item+4001)
                    rad, xi, xierr, DD, DR, RR = np.loadtxt(thisout, unpack=True)
                    plt.scatter(rad, xi*rad**2)
                    plt.axis([-5, 155, 0, 120])
                    plt.xlabel('$\mathrm{\ r \ (Mpc)}$')
                    plt.ylabel(r'$\mathrm{\ r^2 * \xi}$')
                    plt.savefig(os.path.join(drdir, 'qaplots', 'power_spectrum_monopole_{}.png'.format(item+4001)))
                    
            if args.qaplots == '3D_rm':
                for item in range(len(randomslist)):
                    mu, rad, xi, xierr, DD, DR, RR = np.loadtxt(thisout, unpack=True)
                    mono1 = compute_monopole(mu, rad, xi)
                    q1 = compute_quadrupole(mu, rad, xi)
                    hex1 = compute_hexadecapole(mu, rad, xi)
                    # plt.savefig(os.path.join('drdir, 'xi-with-weights.pdf')) 
                    plt.imshow(xi.reshape(50, 40))
                    plotmqh(mono1,q1,hex1,rad)
                    
            if args.qaplots == '3D_ps':
                for item in range(len(randomslist)):
                    pi, sigma, xi, xierr, DD, DR, RR = np.loadtxt(thisout, unpack=True)
                    mono1 = compute_monopole(pi, sigma, xi)
                    q1 = compute_quadrupole(pi, sigma, xi)
                    hex1 = compute_hexadecapole(pi, sigma, xi)
                    # plt.savefig(os.path.join('drdir, 'xi-with-weights.pdf')) 
                    plt.imshow(xi.reshape(50, 40))
                    plotmqh(mono1,q1,hex1,rad)

       
if __name__ == "__main__":
    main()
