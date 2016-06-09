#!/usr/bin/python

from __future__ import division, print_function

import os
import sys
import argparse
import glob
import numpy as np

def main():
    
    # argument parsing
    parser = argparse.ArgumentParser()

    parser.add_argument('--dr', type=str, default='dr11', help='Specify the SDSS data release.')
    parser.add_argument('--type', type=str, default='3D_rm', help='Specify the correlation type used.')
    parser.add_argument('--cov', action='store_true', help='Compute the covariace matrix.')

    args = parser.parse_args()

    # convenience variables
    datadir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'cuteout', args.type)
    outdir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'covariance'

    if args.cov:
    
        mu, rad, xi, xierr, DD, DR, RR = np.loadtxt(thisout, unpack=True)
        covarance = np.cov()
        
if __name__ == "__main__":
    main()
    
