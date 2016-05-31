#!/usr/bin/env python

"""Generate the CUTE input files from the SDSS catalogs and randoms posted at
http://data.sdss3.org/sas/dr11/boss/lss/

"""

# r perp para
# all 600
# monopole vs published
# implement logging
# calculate random weights and reimpliment data weights
# all combinations
# covariance matrix

from __future__ import division, print_function

import os
import sys
import argparse
import logging as log
import numpy as np
from astropy.io import fits
import astropy.cosmology as cosmo


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

def more_cute():
    # for item in randomsdir:
        # randomfile = item
        # do cute
        # write out unique file names
    return end

# def calc_fkp_distance(z):
#     weight = cosmo.comoving_distance(z)


    
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--dr', type=str, default='dr11', help='Specify the SDSS data release.')
    parser.add_argument('--parse', action='store_true', help='Parse the input datafiles.')
    parser.add_argument('--docute', type=str, default=None, help='Run CUTE.')
    parser.add_argument('--qaplots', type=str, default=None, help='Generate QAplots.')

    args = parser.parse_args()

    key = 'LSS_BOSS'
    if key not in os.environ:
        log.fatal('Required ${} environment variable not set'.format(key))
        return 0

    CUTEdir = os.path.join(os.getenv('CUTE'))
    drdir = os.path.join(os.getenv('LSS_BOSS'), args.dr)
    randomsdir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'randoms','*.dat')
    datafile = os.path.join(drdir, args.dr+'_cmass.dat')
    randomfile = os.path.join(drdir, args.dr+'_cmass_random.dat')
    outfile = os.path.join(drdir, 'dr11_2pt_rad.dat')
    paramfile = os.path.join(drdir, 'dr11_rad.param')
    # unique names for output and param files

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

        log.info('Writing {}'.format(datafile))
        np.savetxt(datafile, data)
	
        ra, dec, z, ipoly, wboss, wcp, wzf, veto = \
          np.loadtxt(os.path.join(drdir, 'mock_random_DR11_CMASS_N_PTHALOS_ir4001.dat'), unpack=True)
        keep = np.where(veto==1)[0]
        nobj = len(keep)
	
        rand = np.zeros((nobj,4))
        rand[:,0] = ra[keep]
        rand[:,1] = dec[keep]
        rand[:,2] = z[keep]
        rand[:,3] = wcp[keep]+wzf[keep]-1

        log.info('Writing {}'.format(randomfile))
        np.savetxt(randomfile, rand)

    if args.docute:

        for item in randomsdir:

            randomfile = os.path.join(randomsdir, str(item))

            pfile = open(paramfile,'w')
            pfile.write('data_filename= '+datafile+'\n')
            pfile.write('random_filename= '+randomfile+'\n')
            pfile.write('mask_filename= junk\n')
            pfile.write('z_dist_filename= junk\n')
            pfile.write('output_filename= '+outfile+'\n')
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
                
            if args.docute == '3D_rm':
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
        
            os.system('CUTE '+paramfile)
            
            #os.system('CUTE-noweights '+paramfile)

    if args.qaplots:

        import PIL
        import matplotlib.pyplot as plt
        
        # Make rockin' plots and write out.
        if args.qaplots == 'monopole':
            rad, xi, xierr, DD, DR, RR = np.loadtxt(outfile, unpack=True)
            plt.figure()
            plt.scatter(rad, xi*rad**2)
            plt.axis([-5, 155, 0, 120])
            plt.xlabel('$\mathrm{\ r \ (Mpc)}$')
            plt.ylabel(r'$\mathrm{\ r^2 * \xi}$')
            # plt.savefig(os.path.join('/home/work/projects/lss-boss/dr11', 'xi-with-weights.pdf'))
            plt.show()

        if args.qaplots == '3D_rm':
            mu, rad, xi, xierr, DD, DR, RR = np.loadtxt(outfile, unpack=True)
            mono1 = compute_monopole(mu, rad, xi)
            q1 = compute_quadrupole(mu, rad, xi)
            hex1 = compute_hexadecapole(mu, rad, xi)
            # plt.savefig(os.path.join('/home/work/projects/lss-boss/dr11', 'xi-with-weights.pdf')) 
            plt.imshow(xi.reshape(50, 40))
            plotmqh(mono1,q1,hex1,rad)
            # plt.show()
            # Make 2x2 matrix of images
        
if __name__ == "__main__":
    main()
