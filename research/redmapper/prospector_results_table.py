# script that reads in prospector results and organizes them into a
# table.

import os
import sys
import pdb

import matplotlib.pylab as plt
import numpy as np
from glob import glob
from astropy.table import Table

import h5py
from prospect.io import read_results
import fitsio
import seaborn as sns

from prospector_utilities import logmass2mass
from astropy.table import Table, Column, vstack, hstack

def datadir():
    """Return the top-level data directory."""
    return os.path.join(os.getenv('HOME'), 'stellar-mass')

def read_redmapper():
    """Read the parent redmapper catalog."""
    redfile = os.path.join(os.sep, 'global', 'work', 'projects', 
                    'redmapper', 'redmapper_isedfit_v5.10_centrals.fits.gz')
    print('Reading {}'.format(redfile))
    cat = fitsio.read(redfile, ext=1)
    return cat

prefix = 'redmapper_sdssphot'
these = np.arange(36) # prospector sample objects
cat = read_redmapper()
out = cat[these]

sample = glob('/home/desi2/stellar-mass/redmapper_sdssphot_*_mcmc.h5')

# horrible way of doing it......
mass_results = []
tau_results = []
tage_results = []
zsol_results = []
dust_results = []
for obj in out:
    objprefix = '{0:05}'.format(obj['ISEDFIT_ID'])
    h5file = os.path.join( datadir(), '{}_{}_mcmc.h5'.format(prefix, objprefix) )
    results, guesses, model = read_results.results_from(h5file,model_file=None)

    nwalkers, niter, nparams = results['chain'][:, :, :].shape
    flatchain = results['chain'].reshape(nwalkers * niter, nparams)
    lnp = results['lnprobability'].reshape(nwalkers * niter)
    theta = flatchain[lnp.argmax(), :] # maximum likelihood values

    # make an array with the max liklihood values of each object in the sample. 
    mass_results.append(theta[0])
    tau_results.append(theta[1])
    tage_results.append(theta[2])
    zsol_results.append(theta[3])
    dust_results.append(theta[4])

rownumber = 36 # len(sample)
prospect_table = Table()

prospect_table.add_column(Column(name='objid', length=rownumber))
prospect_table.add_column(Column(name='mass', length=rownumber, dtype='f8'))
prospect_table.add_column(Column(name='tau', length=rownumber, dtype='f8'))
prospect_table.add_column(Column(name='tage', length=rownumber, dtype='f8'))
prospect_table.add_column(Column(name='zsol', length=rownumber, dtype='f8'))
prospect_table.add_column(Column(name='dust', length=rownumber, dtype='f8'))



#prospect_table['objid'] = mass_results
prospect_table['mass'] = mass_results
prospect_table['tau'] = tau_results
prospect_table['tage'] = tage_results
prospect_table['zsol'] = zsol_results
prospect_table['dust'] = dust_results

pdb.set_trace()
print(prospect_table)
