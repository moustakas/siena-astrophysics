#!/usr/bin/env python

"""Wrapper on prospector to derive stellar masses for a subset of the redmapper
sample.

"""
from __future__ import division, print_function

import os
import argparse
from time import time

import numpy as np
import fitsio

def datadir():
    """Return the top-level data directory."""
    return os.path.join(os.getenv('HOME'), 'massivepilot')    

def read_redmapper():
    """Read the parent redmapper catalog."""
    redfile = os.path.join(os.sep, 'global', 'work', 'projects', 
                           'redmapper', 'redmapper_isedfit_v5.10_centrals.fits.gz')
    print('Reading {}'.format(redfile))
    cat = fitsio.read(redfile, ext=1)
    return cat

def read_parent():
    """Read the parent (pilot) massivepilot catalog."""
    fitsfile = os.path.join(datadir(), 'massivepilot-sample.fits')
    print('Reading {}'.format(fitsfile))
    cat = fitsio.read(fitsfile, ext=1)
    return cat

def getobs(cat):
    """Generate the prospector-style "obs" dictionary which contains the input
    photometry, redshift, etc. for a single object.

    cat - fitsio structured array

    """
    from sedpy.observate import load_filters

    sdss = ['sdss_{}0'.format(b) for b in ['u','g','r','i','z']]
    wise = ['wise_{}'.format(b) for b in ['w1','w2']]
    filternames = sdss + wise
    
    mask = (cat['IVARMAGGIES'] > 0) * 1
    obs = {}
    obs['maggies'] = np.squeeze(cat['MAGGIES'])
    with np.errstate(divide='ignore'):
        obs['maggies_unc'] = np.squeeze(1.0/np.sqrt(cat['IVARMAGGIES'])) #[:3, :])
    obs['wavelength'] = None # not fitting spectra
    obs['filternames'] = filternames
    obs['filters'] = load_filters(filternames)
    obs['phot_mask'] = mask  # 1 = good, 0 = bad
    obs['isedfit_id'] = cat['ISEDFIT_ID']
    obs['zred'] = cat['Z']
    
    print('What is "unc"?!?')
    obs['unc'] = obs['maggies_unc']
    obs['spectrum'] = None
    
    return obs

def load_model(zred):
    """Initialize the priors on each free and fixed parameter.

    TBD: Do we need to define priors on dust, fburst, etc., etc.???

    Args:
      zred (float): input galaxy redshift.

    Returns:
      sed (prospect.models.sedmodel.SedModel): SED priors and other stuff.

    """
    from prospect.models import priors, sedmodel
    
    model_params = []
    
    # (Fixed) prior on galaxy redshift.
    model_params.append({'name': 'zred', 'N': 1,
                         'isfree': False,
                         'init': zred,
                         'units': '',
                         'prior_function': priors.tophat,
                         'prior_args': {'mini': 0.0, 'maxi': 4.0}})
    
    # Priors on stellar mass and stellar metallicity.
    model_params.append({'name': 'mass', 'N': 1,
                         'isfree': True,
                         'init': 1e11,
                         'init_disp': 1e10,
                         'units': r'M_\odot',
                         'prior_function': priors.tophat,
                         'prior_args': {'mini': 1e10, 'maxi': 5e12}})

    model_params.append({'name': 'logzsol', 'N': 1,
                         'isfree': False,
                         'init': -0.3,
                         'init_disp': 0.3,
                         'units': r'$\log (Z/Z_\odot)$',
                         'prior_function': priors.tophat,
                         'prior_args': {'mini':-1, 'maxi':0.19}})

    model_params.append({'name': 'mass_units', 'N': 1, # Ben's speed solution.
                        'isfree': False,
                        'init': 'mformed'})
    
    # Priors on SFH type (fixed), tau, and age.
    model_params.append({'name': 'sfh', 'N': 1,
                         'isfree': False,
                         'init':   4, # 4 = delayed tau model
                         'units': 'type'})

    model_params.append({'name': 'tau', 'N': 1,
                         'isfree': True,
                         'init':      1.0,
                         'init_disp': 0.5,
                         'units': 'Gyr',
                         'prior_function': priors.logarithmic,
                         'prior_args': {'mini': 0.1, 'maxi': 5.0}})

    model_params.append({'name': 'tage', 'N': 1,
                        'isfree': True,
                        'init':      10.0,
                        'init_disp':  3.0,
                        'units': 'Gyr',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini': 0.5, 'maxi': 14.0}})
    
    return sedmodel.SedModel(model_params)

def lnprobfn(theta, model, obs, sps, spec_noise=None, phot_noise=None, verbose=True):
    """Define the likelihood function.

    Given a parameter vector and a dictionary of observational data and a model
    object, return the ln of the posterior. This requires that an sps object
    (and if using spectra and gaussian processes, a GP object) be instantiated.

    """
    from prospect.likelihood import lnlike_spec, lnlike_phot, write_log

    lnp_prior = model.prior_product(theta)
    if np.isfinite(lnp_prior):
        
        # Generate mean model
        t1 = time()
        try:
            mu, phot, x = model.mean_model(theta, obs, sps=sps)
        except(ValueError):
            return -np.infty
        d1 = time() - t1

        # Noise modeling
        if spec_noise is not None:
            spec_noise.update(**model.params)
        if phot_noise is not None:
            phot_noise.update(**model.params)
        vectors = {'spec': mu, 'unc': obs['unc'],
                   'sed': model._spec, 'cal': model._speccal,
                   'phot': phot, 'maggies_unc': obs['maggies_unc']}

        # Calculate log-likelihoods
        t2 = time()
        lnp_spec = lnlike_spec(mu, obs=obs, spec_noise=spec_noise, **vectors)
        lnp_phot = lnlike_phot(phot, obs=obs, phot_noise=phot_noise, **vectors)
        d2 = time() - t2
        if verbose:
            write_log(theta, lnp_prior, lnp_spec, lnp_phot, d1, d2)

        return lnp_prior + lnp_phot + lnp_spec

    else:
        return -np.infty

def chisqfn(theta, model, obs, sps):
    """Return the negative of lnprobfn for minimization."""
    return -lnprobfn(theta, model, obs, sps)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--build-sample', action='store_true', help='Build the sample.')
    parser.add_argument('--do-fit', action='store_true', help='Run prospector!')
    parser.add_argument('--qaplots', action='store_true', help='Make some neat plots.')
    parser.add_argument('--verbose', action='store_true', help='Be loquacious.')

    args = parser.parse_args()

    if args.build_sample:

        # Read the parent redmapper catalog, choose a subset of objects and
        # write out.
        cat = read_redmapper()

        # Choose objects with masses from iSEDfit, Kravtsov, and pymorph, but
        # for now just pick three random galaxies.
        these = [3, 4, 5]
        print('Selecting {} galaxies.'.format(len(these)))
        out = cat[these]

        outfile = os.path.join(datadir(), 'massivepilot-sample.fits')
        print('Writing {}'.format(outfile))
        fitsio.write(outfile, out)
        
    if args.do_fit:
        from scipy.optimize import minimize
        
        from prospect.models import model_setup
        from prospect.sources import CSPSpecBasis
        from prospect import fitting

        # Specify the run parameters.
        run_params = {
            'outfile': 'test',
            'debug':   False,
            'nwalkers': 50, # 128,
            'nburn': [32, 32, 64], 
            'niter': 500, # 512,
            'do_powell': True,
            'ftol': 0.5e-5, 'maxfev': 5000,
            'zcontinuous': 1, # interpolate the models in stellar metallicity
            }

        # Specify the minimization method.
        min_method = 'Nelder-Mead'
            
        # Load the default SPS model.
        t0 = time()
        sps = CSPSpecBasis(zcontinuous=run_params['zcontinuous'], compute_vega_mags=False)
        print('Initializing the CSPSpecBasis Class took {:.1f} seconds.'.format(time() - t0))

        # Read the parent sample and loop on each object.
        cat = read_parent()
        
        for obj in cat:
            # Grab the photometry for this object and then initialize the priors
            # and the SED model.
            obs = getobs(obj)
            model = load_model(obs['zred'])

            # Get close to the right answer doing a simple minimization.
            initial_theta = model.rectify_theta(model.initial_theta) # initial parameters
            
            t0 = time()
            min_results = minimize(chisqfn, initial_theta, (model, obs, sps), method=min_method)
            tt = time() - t0

            print('What is edge_trunc!')
            initial_center = fitting.reinitialize(min_results.x, model) # , edge_trunc=run_params.get('edge_trunc', 0.1))
            initial_prob = -1 * min_results['fun']
            
            print('Minimization {} finished in {} seconds'.format(min_method, tt))
            print('best {0} guess: {1}'.format(min_method, initial_theta))
            print('best {0} lnp: {1}'.format(min_method, initial_prob))

            break

            ## reinitialize fit
            #initial_prob = -1 * min_results['fun']
            
    if args.qaplots:
        # Do it
        pass

if __name__ == "__main__":
    main()
