#!/usr/bin/env python

"""Generate the CUTE input files from the SDSS catalogs and randoms posted at
http://data.sdss3.org/sas/dr11/boss/lss/

"""
from __future__ import division, print_function

import os
import sys
import pdb
import logging
import argparse
import glob

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import WMAP7
from matplotlib.colors import LogNorm

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

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

def _dr11_cmass_north_zminmax():
    '''Return minimum and maximum redshifts for this sample.'''
    return 0.43, 0.7

def _parse_speczcat(sample='dr11_cmass_north', clobber=False):
    '''Parse the spectroscopic redshift catalog for a give sample.'''

    sampledir = os.path.join(os.getenv('LSS_CUTE'), sample)

    if sample == 'dr11_cmass_north':
        zmin, zmax = _dr11_cmass_north_zminmax()
        datafile = os.path.join(sampledir, 'cutefiles', '{}_specz.dat'.format(sample))
        
        if not os.path.isfile(datafile) or clobber:

            speczfile = os.path.join(sampledir, 'galaxy_DR11v1_CMASS_North.fits.gz')
            if not os.path.isfile(speczfile):
                log.fatal('Spectroscopic redshift catalog {} not found!'.format(speczfile))
                return 0

            log.info('Reading {}.'.format(speczfile))
            allspecz = fits.getdata(speczfile, 1)
            keep = np.where((allspecz['Z'] > zmin) * (allspecz['Z'] < zmax))[0]
            specz = allspecz[keep]

            log.info('Calculating FKP weights.')
            fkp = calc_fkp_weights(specz['Z'], zmin, zmax)
            data = np.zeros((len(keep), 4))
            data[:, 0] = specz['RA']
            data[:, 1] = specz['DEC']
            data[:, 2] = specz['Z']
            data[:, 3] = fkp*specz['WEIGHT_SYSTOT']*(specz['WEIGHT_NOZ']+specz['WEIGHT_CP']-1)
            
            log.info('Writing {}'.format(datafile))
            np.savetxt(datafile, data)
            
        return datafile
    
    else:
        log.fatal('Unrecognized sample {}.'.format(sample))
        return 0

def _parse_randomcat(sample='dr11_cmass_north', infile=None, outfile=None, clobber=False):
    '''Parse a given random catalog.

    infile/outfile can both be arrays.

    '''
    if infile is None or outfile is None:
        log.fatal('Input and output file names must be provided')
        return 0

    if sample == 'dr11_cmass_north':
        zmin, zmax = _dr11_cmass_north_zminmax()
        for infile1, outfile1 in zip([infile], [outfile]):
            if not os.path.isfile(outfile1) or clobber:
                log.info('Reading {}'.format(infile1))
                ra, dec, z, ipoly, wboss, wcp, wzf, veto = np.loadtxt(infile1, unpack=True)

                keep = np.where(veto == 1)[0]
                rand = np.zeros((len(keep), 4))
                rand[:, 0] = ra[keep]
                rand[:, 1] = dec[keep]
                rand[:, 2] = z[keep]

                #log.info('  Calculating FKP weights.')
                randfkp = calc_fkp_weights(rand[:, 2], zmin, zmax)
                rand[:, 3] = randfkp * (wcp[keep]+wzf[keep]-1)

                log.info('Writing {}'.format(outfile1))
                np.savetxt(outfile1, rand)

    else:
        log.fatal('Unrecognized sample {}.'.format(sample))

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('--sample', type=str, default='dr11_cmass_north', help='Dataset to use (currently only dr11_cmass_north is supported).')
    parser.add_argument('--omegaM', type=float, default='0.3', help='Omega_matter (note: Omega_Lambda = 1-Omega_Matter)')
    parser.add_argument('--w', type=float, default='-1.0', help='w parameter (choose w=-1.0 for cosmological constant)')
    parser.add_argument('--corrtype', type=str, default='monopole', help='Specify correlation type (monopole|3D_ps|3D_rm).')
    parser.add_argument('--nrandom', type=str, default='all', help='Number of random catalogs to use (integer number|all)')
    parser.add_argument('--docute', action='store_true', help='Generate the individual correlation functions using CUTE.')
    parser.add_argument('--qaplots', action='store_true', help='Generate QAplots.')
    parser.add_argument('--clobber', action='store_true', help='Regenerate the parsed data/random files, even if they exist.')

    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    key = 'LSS_CUTE'
    if key not in os.environ:
        log.fatal('Required ${} environment variable not set'.format(key))
        return 0

    # Choose the sample.
    if args.sample == 'dr11_cmass_north':
        sample = args.sample
        zmin, zmax = _dr11_cmass_north_zminmax()
    else:
        log.warning('Unrecognized sample {}'.format(args.sample))
        return 0

    # Establish the top-level directory.
    sampledir = os.path.join(os.getenv('LSS_CUTE'), sample)
    if not os.path.isdir(sampledir):
        log.fatal('Top level directory {} does not exist!'.format(sampledir))
        return 0

    randomdir = os.path.join(sampledir, 'randoms')
    if not os.path.isdir(randomdir):
        log.fatal('Required randoms directory {} does not exist!'.format(randomdir))
        return 0

    # Initialize the cosmological parameters.
    omega_M = args.omegaM
    omega_L = 1-omega_M
    ww = args.w

    # Sensibly hard-code some parameters specific to each correlation function. 
    if args.corrtype == 'monopole':
        log_bin = 0    # logarithmic binning? [should also try 1]
        n_logint = 25  # number of logarithmic bins

        dim1_max = 150  # [Mpc]
        dim1_nbin = 75  # number of bins
        dim2_max = 150  # used? [Mpc]
        dim2_nbin = 75  # number of bins
        
    if args.corrtype == '3D_ps':
        log_bin = 0   # no logarithmic binning
        n_logint = 1  # dummy variable; not used

        dim1_max = 150   # maximum "pi" value [Mpc]
        dim1_nbin = 150  # number of "pi" bins
        dim2_max = 150   # maximum "sigma" value [Mpc]
        dim2_nbin = 150  # number of "sigma" bins

    if args.corrtype == '3D_rm':
        log_bin = 0    # logarithmic binning? [should also try 1]
        n_logint = 25  # number of logarithmic bins

        dim1_max = 150   # maximum "r" value [Mpc]
        dim1_nbin = 150  # number of "r" bins
        dim2_max = 1     # maximum "mu" value [Mpc]
        dim2_nbin = 20   # number of "mu" bins

    ##########
    # Run CUTE.
    if args.docute:

        # Make the subdirectories we need.
        newdir = ['cutefiles']
        for dd in newdir:
            try:
                os.stat(os.path.join(sampledir, dd))
            except:
                os.mkdir(os.path.join(sampledir, dd))
        cutefiledir = os.path.join(sampledir, 'cutefiles')
             
        # Parse the spectroscopic data file (unless it exists).
        speczfile = _parse_speczcat(sample, clobber=args.clobber)

        # Call CUTE using each random catalog in turn.
        allrandomfile = np.array(glob.glob(os.path.join(randomdir, '*')))
        if len(allrandomfile) == 0:
            log.fatal('No random catalogs in {} found!'.format(randomdir))
            return 0

        for ii, randomfile in enumerate(allrandomfile):

            randfile = os.path.join(cutefiledir, '{}_{:05d}.dat'.format(sample, ii+1))
            _parse_randomcat(infile=randomfile, outfile=randfile, clobber=args.clobber)

            # Correlation function output file name.
            outfile = os.path.join(cutefiledir, '{}_{:05d}_{}.dat'.format(sample, ii+1, args.corrtype))

            # Write the CUTE parameter file.
            paramfile = os.path.join(cutefiledir, '{}.param'.format(sample))
            
            pfile = open(paramfile, 'w')
            pfile.write('data_filename= {}\n'.format(speczfile))
            pfile.write('random_filename= {}\n'.format(randfile))
            pfile.write('mask_filename= junk\n')
            pfile.write('z_dist_filename= junk\n')
            pfile.write('output_filename= {}\n'.format(outfile))
            pfile.write('corr_type= {}\n'.format(args.corrtype))
            pfile.write('num_lines= all\n')
            pfile.write('corr_estimator= LS\n')
            pfile.write('input_format= 2\n')
            pfile.write('np_rand_fact= 1\n')
            pfile.write('omega_M= {}\n'.format(omega_M))
            pfile.write('omega_L= {}\n'.format(omega_L))
            pfile.write('w= {}\n'.format(ww))
            pfile.write('radial_aperture= 1\n') # [degrees]
            pfile.write('use_pm= 0\n')
            pfile.write('n_pix_sph= 2048\n')

            pfile.write('log_bin= {}\n'.format(log_bin))
            pfile.write('n_logint= {}\n'.format(n_logint))
            pfile.write('dim1_max= {}\n'.format(dim1_max))
            pfile.write('dim1_nbin= {}\n'.format(dim1_nbin))
            pfile.write('dim2_max= {}\n'.format(dim2_max))
            pfile.write('dim2_nbin= {}\n'.format(dim2_nbin))
            pfile.write('dim3_min= {}\n'.format(zmin))
            pfile.write('dim3_max= {}\n'.format(zmax))
            pfile.write('dim3_nbin= 1\n')
            pfile.close()

            # Run CUTE, passing the newly created parameter file
            os.system('CUTE {}'.format(paramfile))

            pdb.set_trace()

    ##########
    # Run CUTE.
    if args.qaplots:
        
        anderson1 = os.path.join(sampledir, 'Anderson_2013_CMASSDR11_corrfunction_x0x2_prerecon.dat')
        and_rad,and_mono,and_quad = np.loadtxt(anderson1, unpack=True)

        if args.corrtype == 'monopole':
            for item in range(len(randomslist)):
                thisout = outfile+'{}.dat'.format(item+4001)
                rad, xi, xierr, DD, DR, RR = np.loadtxt(thisout, unpack=True)
                plt.scatter(rad, xi*rad**2)
                #plt.axis([-5, 155, 0, 120])
                plt.xlabel('$\mathrm{\ r \ (Mpc)}$')
                plt.ylabel(r'$\mathrm{\ r^2 * \xi}$')
                #plt.savefig(os.path.join(sampledir,'qaplots','power_spectrum_monopole_{}.png'.format(item+4001)))
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
