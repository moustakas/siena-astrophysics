#!/usr/bin/python

from __future__ import division, print_function

import os
import sys
import argparse
import glob

import numpy as np
import matplotlib.pyplot as plt
import pdb

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

    # Read one correlation function and get the 40, 50 dimensions.
    corrfiles = corrfiles[:20] # testing!

    ncorr = len(corrfiles)
    #allcorr = np.zeros((ncorr, 40, 50))
    allcorr = np.zeros((ncorr, 2000))
    for ii, cfile in enumerate(corrfiles):
        print('Reading {}'.format(cfile))
        corr = np.loadtxt(cfile)#.reshape(40, 50, 7)
        allcorr[ii, :] = corr[:,2] # grab xi
        #allcorr[ii, :, :] = corr[:, :, 2] # grab xi

    #cov = np.zeros((40, 50))
    cov = np.zeros((20, 2000))
    meancorr = np.mean(allcorr, axis=1)
    for mm, corrf in enumerate(allcorr):
        for ii in corrf:
            for jj in range(2000):
                cov[mm, jj] = (ii-meancorr[mm])*(corrf[jj]-meancorr[mm])
    pdb.set_trace()
        
    # sys.exit(1)

    # xi = []
    # xibars = []
    # covlist = []
    
    # if args.cov:

    #     for ii in range(len(allfiles)):
    #         xi.append(np.loadtxt(allfiles[ii])[:,2])
    #         xi[ii] = np.reshape(xi[ii], [40,50])
    #         print(ii)
            
    #     for ii in range(len(xi)):
    #         xibars.append(np.mean(xi[ii]))
            
    #     for mm in range(np.shape(xi)[0]): # loop through each file
    #         for ii in range(np.shape(xi)[1]): # loop through each ...
    #             for jj in range(np.shape(xi)[2]): # loop through each ...
    #                 covlist.append((xi[mm][ii]-xibars[mm])*(xi[mm][jj]-xibars[mm]))
                    
if __name__ == "__main__":
    main()
    
