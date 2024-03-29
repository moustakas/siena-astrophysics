#!/usr/bin/python

from __future__ import division, print_function

import os
import sys
import argparse
import glob

import numpy as np
import matplotlib.pyplot as plt

def main():
    
    # argument parsing
    parser = argparse.ArgumentParser()

    parser.add_argument('--dr', type=str, default='dr11', help='Specify the SDSS data release.')
    parser.add_argument('--type', type=str, default='3D_rm', help='Specify the correlation type used.')
    parser.add_argument('--cov', action='store_true', help='Compute the covariace matrix.')

    args = parser.parse_args()

    # convenience variables
    datadir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'cuteout', args.type)
    corrfiles = glob.glob(os.path.join(datadir, '*_fkp_????.dat'))
    outdir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'covariance')

    #corrfiles = corrfiles[:20] # testing!
    thisdimension = len(np.unique(pi))
    thatdimension = len(np.unique(sig))

    ncorr = len(corrfiles)

    xi = np.zeros((ncorr, thisdimension, thatdimension))
    
    for ii, cfile in enumerate(corrfiles):
        print('Reading {}'.format(cfile))
        data = np.loadtxt(cfile)
        xi[ii, :] = data[:,2].reshape(thisdimension,thatdimension) # grab xi and reshape

    rad = np.unique(np.loadtxt(cfile)[:,1])
    mu = np.unique(np.loadtxt(cfile)[:,0])
    cov = np.zeros((40,40))
    xibar = np.mean(xi, axis=0)

    xr = 0.025
    xibar = np.mean(np.mean(xi, axis=2), axis=0)
    mono1 = xr*np.trapz(xi, axis=1)

    for mm in range(ncorr):
        for ii in range(0, 75):
            for jj in range(0, 75):
                # cov[ii,jj] += (xibar[mm, ii] - xii[ii])*(xibar[mm, jj] - xij[jj])
                cov[ii,jj] += (xi[mm, ii] - xibar[ii])*(xi[mm, jj] - xibar[jj])
                
    cov = cov/(ncorr-1)
    aa = np.tile(rad[0:40]*rad[0:40], (1, 40)).reshape(40,40)
                        
if __name__ == "__main__":
    main()
    
