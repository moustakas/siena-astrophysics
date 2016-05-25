#!/usr/bin/env python

"""Generate the CUTE input files from the SDSS catalogs and randoms posted at
http://data.sdss3.org/sas/dr11/boss/lss/

"""

from __future__ import division, print_function

import os
import sys
import argparse
import logging as log

import numpy as np
from astropy.io import fits

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--dr', type=str, default='dr11', help='Specify the SDSS data release.')
    parser.add_argument('--parse', action='store_true', help='Parse the input datafiles.')
    parser.add_argument('--docute', action='store_true', help='Run CUTE.')
    parser.add_argument('--qaplots', action='store_true', help='Generate QAplots.')

    args = parser.parse_args()

    # Set up the logger; basic error checking.
    #log = logging.getLogger()

    key = 'LSS_BOSS'
    if key not in os.environ:
        log.fatal('Required ${} environment variable not set'.format(key))
        return 0
    drdir = os.path.join(os.getenv('LSS_BOSS'), args.dr)

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
	data[:,3] = specz['WEIGHT_FKP']*specz['WEIGHT_SYSTOT']*(specz['WEIGHT_NOZ']+specz['WEIGHT_CP']-1)

        datafile = os.path.join(drdir, dr+'_cmass.dat')
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

        randomfile = os.path.join(drdir, dr+'_cmass_random.dat')
        log.info('Writing {}'.format(randomfile))
	np.savetxt(randomfile, rand)

    if args.docute:
        # Do stuff; write paramfile; call cute us os.system()
        pass
    
    if args.qaplots:
        # Make rockin' plots and write out.
        pass
           
if __name__ == "__main__":
    main()
