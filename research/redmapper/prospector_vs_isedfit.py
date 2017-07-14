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
from astropy.table import Table

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
    #parser.add_argument('--gather-results', action='store_true', help='gather the objects sampled on')
    parser.add_argument('--plottype', type=str, default='mass', help='make a mass vs mass plot.')
    parser.add_argument('--prefix', type=str, default='redmapper_sdssphot', help='String to prepend to I/O files.')
    args = parser.parse_args()

    run_params = {
        'prefix':  args.prefix,
        }

    # if args.gather_results():
        
    # sample = glob(os.path.join('redmapper_sdssphot_*_mcmc.h5'))

    # If gather_results::
    #    Find the set of galaxies fitted using glob
    #    Create an empty Table() with all the columns you want.
    #       mass, logzsol, tage, tau, dust1, filename(?), chi2 (convert from lnp), ...
    #    Loop on each object, read the HDF5 file and populate your table.
    #    Write out the table.
    # If makeplots:
    #    Read previous table.
    #    Read iSEDfit results + match.
    #    Make super awesome plots.
        
    if args.plottype.lower() == 'mass':
        # read in the iSEDfit masses
        these = np.arange(36) # prospector sample objects
        cat = read_redmapper()
        out = cat[these] # match the cat catalog with the prospector sample objects
        # read in the hdf5 files and isedfit catalog values
        pros_results = []
        for obj in out:
            objprefix = '{0:05}'.format(obj['ISEDFIT_ID'])
            h5file = os.path.join( datadir(), '{}_{}_mcmc.h5'.format(run_params['prefix'], objprefix) )
            #h5file = os.path.join( datadir(), '{}_{}_mcmc.h5'.format('test', objprefix) )
            pdb.set_trace()
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
        fig, ax = plt.subplots()
        ax.scatter(log_results, out['MSTAR'], color='b')
        ax.set_xlabel('Prospector log M')
        ax.set_ylabel('iSEDfit log M')
        
        qafile = os.path.join(datadir(), '{}_mass.png'.format(args.prefix))
        fig.savefig(qafile)

        
if __name__ == "__main__":
    main()
