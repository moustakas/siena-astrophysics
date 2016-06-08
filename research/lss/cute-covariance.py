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
    parser.add_argument('--cov', action='store_true', help='Compute the covariace matrix.')

    args = parser.parse_args()

    # convenience variables
    CUTEdir = os.path.join(os.getenv('CUTE'))
    drdir = os.path.join(os.getenv('LSS_BOSS'), args.dr)
    datadir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'cuteout', args.type)

    if args.cov:

        covarance = np.cov()

        
    

if __name__ == "__main__":
    main()
    
