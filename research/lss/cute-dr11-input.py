#!/usr/bin/env python

"""Generate the CUTE input files from the SDSS catalogs and randoms posted at
http://data.sdss3.org/sas/dr11/boss/lss/

"""

from __future__ import division, print_function

import os
import sys

import numpy as np
from astropy.io import fits

dr11dir = os.path.join(os.getenv('IM_ARCHIVE_DIR'), 'projects', 'boss-lss')

allspecz = fits.getdata(os.path.join(dr11dir, 'galaxy_DR11v1_CMASS_North.fits.gz'), 1)
keep = np.where((allspecz['Z']>0.43)*(allspecz['Z']<0.7))[0]

specz = allspecz[keep]
ngal = len(keep)

data = np.zeros((ngal,4))
data[:,0] = specz['RA']
data[:,1] = specz['DEC']
data[:,2] = specz['Z']
data[:,3] = specz['WEIGHT_FKP']*specz['WEIGHT_SYSTOT']*(specz['WEIGHT_NOZ']+specz['WEIGHT_CP']-1)

np.savetxt(os.path.join(dr11dir, 'dr11_cmass.dat'), data)


ra, dec, z, ipoly, wboss, wcp, wzf, veto = \
  np.loadtxt(os.path.join(dr11dir, 'mock_random_DR11_CMASS_N_PTHALOS_ir4001.dat'), unpack=True)
keep = np.where(veto==1)[0]
nobj = len(keep)

rand = np.zeros((nobj,4))
rand[:,0] = ra[keep]
rand[:,1] = dec[keep]
rand[:,2] = z[keep]
rand[:,3] = wcp[keep]+wzf[keep]-1

np.savetxt(os.path.join(dr11dir, 'dr11_cmass_random.dat'), rand)
