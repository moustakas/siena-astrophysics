#!/usr/bin/python

from __future__ import division, print_function

import os
import sys
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
import pylab

def main():
    
    # argument parsing
    parser = argparse.ArgumentParser()

    parser.add_argument('--dr', type=str, default='dr11', help='Specify the SDSS data release.')
    parser.add_argument('--type', type=str, default='3D_rm', help='Specify the correlation type used.')
    parser.add_argument('--cov', action='store_true', help='Compute the covariace matrix.')

    args = parser.parse_args()

    # convenience variables
    datadir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'cuteout', args.type)
    allfiles = glob.glob(os.path.join(datadir, '*.dat'))
    outdir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'covariance')

    xi = []
    xibars = []
    covlist = []
    
    if args.cov:

        for ii in range(len(allfiles)):
            xi.append(np.loadtxt(allfiles[ii])[:,2])
            xi[ii] = np.reshape(xi[ii], [40,50])
            print(ii)

        for ii in range(len(xi)):
            xibars.append(np.mean(xi[ii]))

        for mm in range(len(allfiles)):
            for ii in range(np.shape(xi)[1]):
                for jj in range(np.shape(xi)[2]):
                    covlist.append((xi[mm][ii]-xibars[mm])*(xi[mm][jj]-xibars[mm]))

if __name__ == "__main__":
    main()
    
