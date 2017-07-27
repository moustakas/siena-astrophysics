"""Script that reads in prospector results and organizes them into a table.

"""
import os, sys, pdb
import argparse

import numpy as np
from astropy.table import Table, Column

def main():

    from glob import glob
    from prospect.io import read_results
    from prospector_utilities import datadir, logmass2mass

    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', type=str, default='redmapper_sdssphot', help='Prefix for I/O files.')
    args = parser.parse_args()

    datadir = datadir()

    # Gather the list of fitted objects.
    h5files = np.array( sorted(glob( os.path.join(datadir, '{}_*_mcmc.h5'.format(args.prefix)))) )
    nobj = len(h5files)
    if nobj == 0:
        print('No files found!')
        sys.exit(1)

    # Initialize the output table.
    prospect_table = Table()
    prospect_table.add_column(Column(name='isedfit_id', length=nobj, dtype='i8'))
    prospect_table.add_column(Column(name='logmass', length=nobj, dtype='f4'))
    prospect_table.add_column(Column(name='logzsol', length=nobj, dtype='f4'))
    prospect_table.add_column(Column(name='dust2', length=nobj, dtype='f4'))
    prospect_table.add_column(Column(name='tau', length=nobj, dtype='f4'))
    prospect_table.add_column(Column(name='tage', length=nobj, dtype='f4'))
    prospect_table.add_column(Column(name='chi2', length=nobj, dtype='f4'))

    for ii, h5file in enumerate(h5files):
        print('Reading results {}/{}: {}'.format(ii+1, nobj, h5file))
        results, guesses, model = read_results.results_from(h5file)

        # Get the maximum likelihood fitting results.
        nwalkers, niter, nparams = results['chain'][:, :, :].shape
        flatchain = results['chain'].reshape(nwalkers * niter, nparams)
        lnp = results['lnprobability'].reshape(nwalkers * niter)
        mlindx = lnp.argmax()
        theta = flatchain[mlindx, :]

        isedfit_id = results['obs']['isedfit_id']

        prospect_table['isedfit_id'][ii] = isedfit_id
        prospect_table['logmass'][ii] = theta[0]
        prospect_table['logzsol'][ii] = theta[1]
        prospect_table['dust2'][ii] = theta[2]
        prospect_table['tau'][ii] = theta[3]
        prospect_table['tage'][ii] = theta[4]

        prob = np.exp(lnp - lnp.max())
        prob /= prob.sum()

        prospect_table['chi2'][ii] = - np.log(prob[mlindx])

    print(prospect_table)

    outfile = os.path.join(datadir, '{}_prospector.fits'.format(args.prefix))
    print('Writing {}'.format(outfile))
    prospect_table.write(outfile, overwrite=True)

    ## Read the parent isedfit catalog and pull out the fitted objects.  Assume
    ## everything was fitted and is properly sorted.
    #parent, _ = read_parent(args.prefix, datadir=datadir)
    #
    #outfile = os.path.join(datadir, '{}_isedfit.fits'.format(args.prefix))
    #print('Writing {}'.format(outfile))
    #Table(parent).write(outfile, overwrite=True)

if __name__ == '__main__':
    main()
