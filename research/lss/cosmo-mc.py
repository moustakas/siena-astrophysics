#!/usr/bin/python

from __future__ import division, print_function

import os
import sys
import argparse
import glob
import logging as log
import numpy as np

def main():

    # create directory and environment variable
    CosmoMC = os.path.join(os.getenv('COSMOMC'))
    outdir = os.path.join(Cosmo, 'cosmoMCout')
    paramfile = os.path.join(Cosmo, 'param')

    parser = argparse.ArgumentParser()

    parser.add_argument('--docosmo', action='store_true', help='Run CosmoMC.')

    args = parser.parse_args()

    key = 'COSMOMC'
    if key not in os.environ:
        log.fatal('Required ${} environment variable not set'.format(key))
        return 0

    if args.docosmo:
        # write the parameter file
        pfile = open(paramfile, 'w')
        pfile.write('\n')
        pfile.close()            

if __name__ == "__main__":
    main()
