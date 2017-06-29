'''
Script to read prospector/emcee minimization results and compare them with Dr. Moustakas' iSEDfit results
'''

import os
import sys
import pdb
import argparse

import matplotlib as plt
import numpy as np
from glob import glob

import h5py
from prospect.io import read_results
import fitsio

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


def main():
    # gather information from hdf5 files and make a plot of results vs isedfit results
    parser = argparse.ArgumentParser()
    parser.add_argument('--plottype', action='store_true', help='choose axis for plot')
    parser.add_argument('--prefix', type=str, default='redmapper_sdssphot', help='String to prepend to I/O files.')
    args = parser.parse_args()

    run_params = {
        'prefix':  args.prefix,
        }
        
    if args.plottype:
        these = np.arange(4)
        cat = read_redmapper()
        out = cat[these]
        # read in the hdf5 files and isedfit catalog values
        pros_results = []
        for obj in cat:
            objprefix = '{0:05}'.format(obj['ISEDFIT_ID'])
            h5file = os.path.join( datadir(), '{}_{}_mcmc.h5'.format(run_params['prefix'], objprefix) )
            #h5file = os.path.join( datadir(), '{}_{}_mcmc.h5'.format('test', objprefix) )
            results, guesses, model = read_results.results_from(h5file,model_file=None)
        
        
            nwalkers, niter, nparams = results['chain'][:, :, :].shape
            flatchain = results['chain'].reshape(nwalkers * niter, nparams)
            lnp = results['lnprobability'].reshape(nwalkers * niter)
            theta = flatchain[lnp.argmax(), :] # maximum likelihood values

            
            pros_results.append(theta[0])
        # make plots of our results vs isedfit. Default is mass vs mass.
        plt.plot(np.log(pros_results), cat['MSTAR'])
        qafile = os.path.join(datadir(), '{}_{}_compare.png'.format(args.prefix, objprefix))
        plt.savefig(qafile)

        
if __name__ == "__main__":
    main()
