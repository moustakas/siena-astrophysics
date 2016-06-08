#!/usr/bin/python

from __future__ import division, print_function

import os
import sys
import argparse
import glob
import numpy as np

def main():

    CUTEdir = os.path.join(os.getenv('CUTE'))
    drdir = os.path.join(os.getenv('LSS_BOSS'), args.dr)
    randomsdir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'randoms')
    datafile = os.path.join(drdir, args.dr+'_cmass.dat')
    randomfile = os.path.join(drdir, 'parsed', args.dr+'_cmass_random')
    outfile = os.path.join(drdir, 'cuteout', 'dr11_2pt_'+args.type+'_')
    paramfile = os.path.join(drdir, 'param', 'dr11_'+args.type+'_')
    randomslist = glob.glob(os.path.join(randomsdir, '*.dat'))

if __name__ == "__main__":
    main()
    
