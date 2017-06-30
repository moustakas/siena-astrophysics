'''
Script to read prospector/emcee minimization results and compare them with Dr. Moustakas' iSEDfit results
'''

import os
import sys
import pdb
import argparse

import matplotlib.pylab as plt
import numpy as np
from glob import glob

import h5py
from prospect.io import read_results
import fitsio
import seaborn as sns

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
        # read in the iSEDfit masses
        these = np.arange(4) # prospector sample objects
        cat = read_redmapper()
        out = cat[these] # match the cat catalog with the prospector sample objects
        # read in the hdf5 files and isedfit catalog values
        pros_results = []
        for obj in out:
            objprefix = '{0:05}'.format(obj['ISEDFIT_ID'])
            h5file = os.path.join( datadir(), '{}_{}_mcmc.h5'.format(run_params['prefix'], objprefix) )
            #h5file = os.path.join( datadir(), '{}_{}_mcmc.h5'.format('test', objprefix) )
            results, guesses, model = read_results.results_from(h5file,model_file=None) #just care about the results
        
            # calculate the max likehood value
            nwalkers, niter, nparams = results['chain'][:, :, :].shape
            flatchain = results['chain'].reshape(nwalkers * niter, nparams)
            lnp = results['lnprobability'].reshape(nwalkers * niter)
            theta = flatchain[lnp.argmax(), :] # maximum likelihood values

            # make an array with the max liklihood values of each object in the sample. 
            pros_results.append(theta[0])
            
        # make plots of our results vs isedfit. Default is mass vs mass.
        pros_results = np.asarray(pros_results)
        log_results = np.log10(pros_results)

        #seaborn.regplot(log_results, out['MSTAR'], x = 'prospector (logm)', y = 'iSEDfit (logm)') 
        plt.figure()
        plt.plot(log_results, out['MSTAR'], color='b')
        plt.xlabel('Prospector log M')
        plt.ylabel('iSEDfit log M')
        plt.title('Prospector vs. iSEDfit')
        
        qafile = os.path.join(datadir(), '{}_compare.png'.format(args.prefix))
        plt.savefig(qafile)

        
if __name__ == "__main__":
    main()
