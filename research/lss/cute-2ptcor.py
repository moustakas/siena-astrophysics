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
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.cosmology import WMAP7
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

logging.basicConfig(level=logging.INFO)
log = logging.getLogger()

def plotmqh(monopole,quadrupole,hexadecapole,rrange):
    plt.plot(rrange, monopole*rrange**2, 'ko')

def compute_monopole(mu, r, xirm):
    '''Compute the monopole from the xi(r, mu) correlation function.'''
    xirm = xirm*1.0
    Bxirm = np.reshape(xirm, [150, 20]) # generalize
    #xr = 0.025 # find out about this factor (factor of 1/40)
    xr = 0.5
    monopole = xr*np.trapz(Bxirm, axis=1) # go through the math of why
    return monopole

def compute_quadrupole(mu, r, xirm):
    '''Compute the quadrapole from the xi(r, mu) correlation function.'''
    xirm = xirm*(3*(mu*mu-1.0)*(5.0/2.0))
    Bxirm = np.reshape(xirm, [150, 20])
    xr = 0.025
    quadrupole = xr*np.trapz(Bxirm)
    return quadrupole

def compute_hexadecapole(mu, r, xirm):
    '''Compute the hexadecapole from the xi(r, mu) correlation function.'''
    xirm = xirm*(35*(mu*mu*mu*mu)-(30*mu*mu)+3)/8
    Bxirm = np.reshape(xirm,[40,50])
    xr = 0.025
    hexadecapole = xr*np.trapz(Bxirm)
    return hexadecapole

def calc_fkp_weights(z, zmin, zmax, area=1.0):
    '''Compute the FKP statistical weights.'''

    NRB = 200 
    dz = zmax - zmin

    omega = area/(4*np.pi*(180.0/np.pi)**2)
    
    red_interval = dz/NRB
    red_markers = []
    red_vol = []
    bin_num = []
    bin_sum = []
    wfkp = []

    for ii in range(NRB + 1):
        red_markers.append(zmin + ii*red_interval) 
        if ii >= 1:
            red_vol.append((4/3)*np.pi*((WMAP7.comoving_distance(red_markers[ii]).value*0.704)**3-
                                        (WMAP7.comoving_distance(red_markers[ii-1]).value*0.704)**3)
                                        *omega)
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
    '''Return minimum and maximum redshifts and the survey area (deg^2).'''
    return 0.43, 0.7, 6308.0

def _literature(author='anderson', sample='dr11_cmass_north'):
    '''Wrapper function to read various published correlation functions.'''

    sampledir = os.path.join(os.getenv('LSS_CUTE'), sample)

    # Anderson+13 monopole and quadrapole
    if author.lower() == 'anderson':
        litfile = os.path.join(sampledir, 'Anderson_2013_CMASSDR11_corrfunction_x0x2_prerecon.dat')
        if os.path.isfile(litfile):
            rad, mono, quad = np.loadtxt(litfile, unpack=True)
            return rad, mono, quad
        else:
            log.warning('Published correlation function {} not found!'.format(litfile))

def _parse_speczcat(sample='dr11_cmass_north', clobber=False):
    '''Parse the spectroscopic redshift catalog for a give sample.'''

    sampledir = os.path.join(os.getenv('LSS_CUTE'), sample)

    if sample == 'dr11_cmass_north':
        zmin, zmax, area = _dr11_cmass_north_zminmax()
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
            fkp = calc_fkp_weights(specz['Z'], zmin, zmax, area)
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
        zmin, zmax, area = _dr11_cmass_north_zminmax()
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
                randfkp = calc_fkp_weights(rand[:, 2], zmin, zmax, area)
                rand[:, 3] = randfkp * (wcp[keep]+wzf[keep]-1)

                log.info('Writing {}'.format(outfile1))
                np.savetxt(outfile1, rand)
                return 1

    else:
        log.fatal('Unrecognized sample {}.'.format(sample))
        return 0

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
        zmin, zmax, area = _dr11_cmass_north_zminmax()
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
        
    elif args.corrtype == '3D_ps':
        log_bin = 0   # no logarithmic binning
        n_logint = 1  # dummy variable; not used

        dim1_max = 150   # maximum "pi" value [Mpc]
        dim1_nbin = 150  # number of "pi" bins
        dim2_max = 150   # maximum "sigma" value [Mpc]
        dim2_nbin = 150  # number of "sigma" bins
        
    elif args.corrtype == '3D_rm':
        log_bin = 0    # logarithmic binning? [should also try 1]
        n_logint = 25  # number of logarithmic bins

        dim1_max = 150   # maximum "r" value [Mpc]
        dim1_nbin = 150  # number of "r" bins
        dim2_max = 1     # maximum "mu" value [Mpc]
        dim2_nbin = 20   # number of "mu" bins
    else:
        log.fatal('Unrecognized or unsupported correlation type {}'.format(args.corrtype))
        return 0

    # Make the subdirectories we need.
    newdir = ['cutefiles', 'qaplots']
    for dd in newdir:
        try:
            os.stat(os.path.join(sampledir, dd))
        except:
            os.mkdir(os.path.join(sampledir, dd))
    qadir = os.path.join(sampledir, 'qaplots')
    cutefiledir = os.path.join(sampledir, 'cutefiles')
             
    ##########
    # Run CUTE.
    if args.docute:

        # Parse the spectroscopic data file (unless it exists).
        speczfile = _parse_speczcat(sample, clobber=args.clobber)
        if type(speczfile) == int:
            return

        # Call CUTE using each random catalog in turn, optionally restricting to
        # a smaller number of random catalogs.
        allrandomfile = glob(os.path.join(randomdir, '*'))
        if len(allrandomfile) == 0:
            log.fatal('No random catalogs in {} found!'.format(randomdir))
            return 0
        if args.nrandom != 'all':
            allrandomfile = allrandomfile[:int(args.nrandom)]

        for ii, randomfile in enumerate(allrandomfile):

            randfile = os.path.join(cutefiledir, '{}_{:05d}.dat'.format(sample, ii+1))
            check = _parse_randomcat(infile=randomfile, outfile=randfile, clobber=args.clobber)
            if check == 0:
                return
            
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

    ##########
    # Generate QAplots.

    if args.qaplots:
        log.info('Building {} QAplots.'.format(args.corrtype))

        allcorrfile = glob(os.path.join(cutefiledir, '{}_?????_{}.dat'.format(sample, args.corrtype)))
        ncorr = len(allcorrfile)

        if args.corrtype == 'monopole':
            allxi = np.zeros((ncorr, dim1_nbin))
            for ii, corrfile in enumerate(allcorrfile):
                if ((ii + 1) % 10) == 0:
                    log.info('Reading correlation function {}/{}'.format(ii+1, ncorr))
                rad, xi, xierr, DD, DR, RR = np.loadtxt(corrfile, unpack=True)
                allxi[ii, :] = xi

            # Compare with Anderson+
            andrad, andmono, andquad = _literature(author='anderson', sample=sample)
                
            xibar = np.mean(allxi, axis=0)

            qafile = os.path.join(qadir, '{}_monopole.pdf'.format(sample))
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(rad, xibar*rad*rad, label='Siena Average Monopole')
            ax.scatter(andrad, andmono*andrad*andrad, marker='s', color='orange', label='Anderson+13 Monopole')
            #ax.scatter(andrad, -andquad*andrad*andrad, marker='s', color='red', label='Anderson+13 Quadrapole')
            ax.set_xlabel(r'$r$ (Mpc / h)')
            ax.set_ylabel(r'$r^2 \xi$')
            ax.legend(loc='upper right', frameon=None)
            ax.margins(0.05)

            plt.subplots_adjust(bottom=0.15, top=0.88)

            log.info('Writing {}'.format(qafile))
            plt.savefig(qafile)

        if args.corrtype == '3D_ps':
            allxi = np.zeros((ncorr, dim1_nbin, dim2_nbin))
            for ii, corrfile in enumerate(allcorrfile):
                if ((ii + 1) % 10) == 0:
                    log.info('Reading correlation function {}/{}'.format(ii+1, ncorr))
                pi, sigma, xi, xierr, DD, DR, RR = np.loadtxt(corrfile, unpack=True)
                allxi[ii, :, :] = xi.reshape(dim1_nbin, dim2_nbin)

            xibar = np.mean(allxi, axis=0)

            bigxi = np.zeros((dim1_nbin*2, dim2_nbin*2))
            bigxi[:dim1_nbin, :dim2_nbin] = np.fliplr(np.flipud(xibar))
            bigxi[dim1_nbin:2*dim1_nbin, :dim2_nbin] = np.fliplr(xibar)
            bigxi[:dim1_nbin, dim2_nbin:2*dim2_nbin] = np.flipud(xibar)
            bigxi[dim1_nbin:2*dim1_nbin, dim2_nbin:2*dim2_nbin] = xibar

            # Fragile!
            pi2d = np.concatenate((np.arange(-dim2_nbin, 0, 1), np.arange(1, dim2_nbin+1, 1)))
            pi2d = np.tile(pi2d, (1, 2*dim2_nbin)).reshape(2*dim2_nbin, 2*dim2_nbin)
            sig2d = np.rot90(pi2d)
            
            qafile = os.path.join(qadir, '{}_2d_ps.pdf'.format(sample))
            fig, ax = plt.subplots(figsize=(8, 8))
            im = plt.pcolor(sig2d, pi2d, bigxi, norm=LogNorm(), cmap='Blues_r')
            #C = plt.contour(sig2d, pi2d, bigxi, colors='k')
            
            ax.set_aspect('equal')
            ax.set_xlabel(r'$\sigma$ (Mpc / h)', fontsize=16)
            ax.set_ylabel(r'$\pi$ (Mpc / h)', fontsize=16)
            ax.set_xlim((-130, 130))
            ax.set_ylim((-130, 130))
            ax.margins(0)

            div = make_axes_locatable(ax)
            cax = div.append_axes('right', size='5%', pad=0.2)
            cbar = fig.colorbar(im, cax=cax, format='%.4g')
            cbar.ax.set_ylabel(r'$\xi(\pi, \sigma)$', fontsize=16)

            plt.subplots_adjust(right=0.85)

            log.info('Writing {}'.format(qafile))
            plt.savefig(qafile)

            pdb.set_trace()

        if args.corrtype == '3D_rm':
            allmono = np.zeros((ncorr, dim1_nbin))
            allquad = np.zeros((ncorr, dim1_nbin))
            for ii, corrfile in enumerate(allcorrfile[:10]):
                if ((ii + 1) % 10) == 0:
                    log.info('Reading correlation function {}/{}'.format(ii+1, ncorr))
                mu, rad, xi, xierr, DD, DR, RR = np.loadtxt(corrfile, unpack=True)
                allmono[ii, :] = compute_monopole(mu, rad, xi)
                #allquad[ii, :] = compute_quadrapole(mu, rad, xi)

            rad = np.unique(rad)
            monobar = np.mean(allmono, axis=0)
            #quadbar = np.mean(allquad, axis=0)
            #pdb.set_trace()

            # Compare with Anderson+
            andrad, andmono, andquad = _literature(author='anderson', sample=sample)

            qafile = os.path.join(qadir, '{}_rm_monopole.pdf'.format(sample))
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(rad, monobar*rad*rad, label='Siena Average Monopole')
            ax.scatter(andrad, andmono*andrad*andrad, marker='s', color='orange', label='Anderson+13 Monopole')
            #ax.scatter(andrad, -andquad*andrad*andrad, marker='s', color='red', label='Anderson+13 Quadrapole')
            ax.set_xlabel(r'$r$ (Mpc / h)')
            ax.set_ylabel(r'$r^2 \xi$')
            ax.legend(loc='upper right', frameon=None)
            ax.margins(0.05)

            plt.subplots_adjust(bottom=0.15, top=0.88)

            log.info('Writing {}'.format(qafile))
            plt.savefig(qafile)

if __name__ == "__main__":
    main()
