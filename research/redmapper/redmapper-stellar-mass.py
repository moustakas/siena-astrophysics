#!/usr/bin/env python

"""Wrapper on prospector to derive stellar masses for a subset of the redmapper
sample.

"""
from __future__ import division, print_function

import os
import sys
import pdb
import argparse
import pickle
from time import time

import numpy as np
import fitsio
import matplotlib.pylab as plt

from prospect.sources import CSPSpecBasis

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

def read_parent():
    """Read the parent (pilot) catalog."""
    fitsfile = os.path.join(datadir(), 'pilot-sample.fits')
    print('Reading {}'.format(fitsfile))
    cat = fitsio.read(fitsfile, ext=1)
    return cat

def getobs(cat):
    """Generate the prospector-style "obs" dictionary which contains the input
    photometry, redshift, etc. for a single object.

    cat - fitsio structured array

    """
    from sedpy.observate import load_filters

    obs = {} 

    # Photometric bandpasses
    sdss = ['sdss_{}0'.format(b) for b in ['u','g','r','i','z']]
    wise = ['wise_{}'.format(b) for b in ['w1','w2']]
    filternames = sdss + wise
    
    obs['filternames'] = filternames
    obs['filters'] = load_filters(filternames)

    # Input photometry
    obs['maggies'] = np.squeeze(cat['MAGGIES'])
    mask = (cat['IVARMAGGIES'] > 0) * 1
    with np.errstate(divide='ignore'):
        obs['maggies_unc'] = np.squeeze(1.0/np.sqrt(cat['IVARMAGGIES'])) #[:3, :])
    obs['phot_mask'] = mask  # 1 = good, 0 = bad

    # Input spectroscopy (none for this dataset)
    obs['wavelength'] = None
    obs['spectrum'] = None
    obs['unc'] = None

    # Use initial values based on the iSEDfit results.
    obs['zred'] = cat['Z'] # redshift
    obs['mass'] = 10**cat['MSTAR'] # stellar mass
    obs['logzsol'] = np.log10(cat['ZMETAL']) # stellar metallicity
    obs['tage'] = cat['AGE']  # age
    obs['tau'] = cat['TAU']   # tau (for a delayed SFH)
    obs['dust2'] = 0.05

    # Additional informational keys.
    obs['isedfit_id'] = cat['ISEDFIT_ID']
    
    return obs

def load_model(zred=0.0, mass=1e11, logzsol=0.0, tage=12.0, tau=1.0, dust2=0.1):
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
                         'init':      mass, 
                         'init_disp': 5e11, 
                         'units': r'M_\odot',
                         'prior_function': priors.tophat,
                         'prior_args': {'mini': 1e9, 'maxi': 1e13}})

    model_params.append({'name': 'logzsol', 'N': 1,
                         'isfree': False,
                         'init': logzsol,
                         'init_disp': 0.3, # dex
                         'units': r'$\log (Z/Z_\odot)$',
                         'prior_function': priors.tophat,
                         'prior_args': {'mini': -1.0, 'maxi': 0.19}})

    model_params.append({'name': 'mass_units', 'N': 1, # Ben's speed solution.
                        'isfree': False,
                        'init': 'mformed'})

    # Priors on dust.
    model_params.append({'name': 'dust2', 'N': 1,
                        'isfree': False,
                        'init':     dust2,
                        'reinit':   True,
                        'init_disp': 0.3,
                        'units': '',
                        'prior_function': priors.tophat,
                        'prior_args': {'mini': 0.0, 'maxi': 2.0}})
    
    # Priors on SFH type (fixed), tau, and age.
    model_params.append({'name': 'sfh', 'N': 1,
                         'isfree': False,
                         'init':   4, # 4 = delayed tau model
                         'units': 'type'})

    model_params.append({'name': 'tau', 'N': 1,
                         'isfree': True,
                         'init':      tau,
                         'init_disp': 10.0,
                         'units': 'Gyr',
                         'prior_function': priors.logarithmic,
                         'prior_args': {'mini': 0.1, 'maxi': 5.0}})
                         #'prior_function': priors.tophat,
                         #'prior_args': {'mini': 0.01, 'maxi': 10.0}})

    model_params.append( {
        'name':   'tage',
        'N':       1,
        'isfree':    True,
        'init':      tage,
        'init_disp':  3.0,
        'units':       'Gyr',
        'prior_function': priors.tophat,
        'prior_args': {'mini': 0.5, 'maxi': 15.0}
        } )
    
    return sedmodel.SedModel(model_params)

def lnprobfn(theta, model, obs, sps, spec_noise=None, phot_noise=None, verbose=False):
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
    parser.add_argument('--prefix', type=str, default='test', help='String to prepend to I/O files.')
    parser.add_argument('--build-sample', action='store_true', help='Build the sample.')
    parser.add_argument('--min-method', default='Powell', type=str,
                        help='Method to use for initial minimization.')
    parser.add_argument('--do-fit', action='store_true', help='Run prospector!')
    parser.add_argument('--qaplots', action='store_true', help='Make some neat plots.')
    parser.add_argument('--threads', default=16, help='Number of cores to use concurrently.')
    parser.add_argument('--remake-sps', action='store_true', help="""Remake (and write out) the
                           CSPSpecBasis object, otherwise read it from disk.""")
    parser.add_argument('--verbose', action='store_true', help='Be loquacious.')

    args = parser.parse_args()

    if args.build_sample:

        # Read the parent redmapper catalog, choose a subset of objects and
        # write out.
        cat = read_redmapper()

        # Choose objects with masses from iSEDfit, Kravtsov, and pymorph, but
        # for now just pick three random galaxies.
        these = [500, 600, 700]
        print('Selecting {} galaxies.'.format(len(these)))
        out = cat[these]

        outfile = os.path.join(datadir(), 'pilot-sample.fits')
        print('Writing {}'.format(outfile))
        fitsio.write(outfile, out, clobber=True)

    if args.do_fit:
        import h5py
        import emcee
        from scipy.optimize import minimize
        
        from prospect.models import model_setup
        from prospect import fitting
        from prospect.io import write_results

        # Specify the run parameters.
        run_params = {
            'prefix': args.prefix,
            'debug':   False,
            'nwalkers': 16, # 128,
            'nburn': [32, 32, 64], 
            'niter': 128, # 512,
            # set the powell convergence criteria 
            'do_powell': False,
            #'ftol': 0.5e-5, 'maxfev': 6000,
            'compute_vega_mags': False,
            'vactoair_flag':      True, # use wavelengths in air
            'zcontinuous': 1,           # interpolate in metallicity
            }

        # Load or generate the default SPS object.
        spsfile = os.path.join( datadir(), '{}_CSPSpecBasis.pickle'.format(args.prefix) )
        if os.path.isfile(spsfile) and not args.remake_sps:
            print('Reading pickled CSPSpecBasis object from {}.'.format(spsfile))
            sps = pickle.load(open(spsfile, 'rb'))
        else:
            t0 = time()
            print('Initializing the CSPSpecBasis object...')
            sps = CSPSpecBasis(zcontinuous=run_params['zcontinuous'],
                               compute_vega_mags=run_params['compute_vega_mags'])
            print('...took {:.1f} seconds.'.format(time() - t0))
            sps.params['vactoair_flag'] = run_params['vactoair_flag']
            sps.update()
            
            print('Pickling CSPSpecBasis to file {}.'.format(spsfile))
            pickle.dump(sps, open(spsfile, 'wb'))

        # Read the parent sample and loop on each object.
        cat = read_parent()
        for obj in cat:
            objprefix = '{0:05}'.format(obj['ISEDFIT_ID'])
            
            # Grab the photometry for this object and then initialize the priors
            # and the SED model.
            obs = getobs(obj)
            model = load_model(zred=obs['zred'], mass=obs['mass'], logzsol=obs['logzsol'],
                               tage=obs['tage'], tau=obs['tau'], dust2=obs['dust2'])

            # Get close to the right answer doing a simple minimization.
            initial_theta = model.rectify_theta(model.initial_theta)

            if False:
                t0 = time()
                min_results = minimize(chisqfn, initial_theta, (model, obs, sps),
                                       method=args.min_method)
                pdur = time() - t0

                print('What is edge_trunc?')
                initial_center = fitting.reinitialize(
                    min_results.x, model, edge_trunc=run_params.get('edge_trunc', 0.1)
                    )
                initial_prob = -1 * min_results['fun']
                
                print('Minimization {} finished in {} seconds'.format(args.min_method, pdur))
                print('best {0} guess: {1}'.format(args.min_method, initial_center))
                print('best {0} lnp: {1}'.format(args.min_method, initial_prob))
            else:
                print('Skipping initial Powell minimization because it borks.')
                min_results, pdur = None, None
                initial_center = initial_theta.copy()
                initial_prob = lnprobfn(initial_center, model, obs, sps)

            hfilename = os.path.join( datadir(), '{}_{}_mcmc.h5'.format(
                run_params['prefix'], objprefix) )
            if os.path.isfile(hfilename):
                os.remove(hfilename)
            
            hfile = h5py.File(hfilename, 'a')
            print('Writing to file {}'.format(hfilename))

            write_results.write_h5_header(hfile, run_params, model)
            write_results.write_obs_to_h5(hfile, obs)
            
            fout = sys.stdout
            fnull = open(os.devnull, 'w')
            sys.stdout = fnull

            tstart = time()
            out = fitting.run_emcee_sampler(lnprobfn, initial_center, model, threads=args.threads, 
                                            initial_prob=initial_prob, hdf5=hfile, nwalkers=run_params.get('nwalkers'),
                                            nburn=run_params.get('nburn'), niter=run_params.get('niter'),
                                            postargs=(model, obs, sps))
            esampler, burn_p0, burn_prob0 = out
            edur = time() - tstart

            sys.stdout = fout
            print('done emcee in {}s'.format(edur))
            
            # Write out more.
        
            write_results.write_pickles(run_params, model, obs, esampler, min_results,
                                        outroot='{}_{}'.format(run_params['prefix'], objprefix),
                                        toptimize=pdur, tsample=edur,
                                        sampling_initial_center=initial_center,
                                        post_burnin_center=burn_p0,
                                        post_burnin_prob=burn_prob0)
            write_results.write_hdf5(hfilename, run_params, model, obs, esampler, 
                                     min_results, toptimize=pdur, tsample=edur,
                                     sampling_initial_center=initial_center,
                                     post_burnin_center=burn_p0,
                                     post_burnin_prob=burn_prob0)
            

    if args.qaplots:        
        import h5py
        from prospect.io import read_results
        from prospector_plot_utilities import param_evol, subtriangle

        import seaborn as sns
        sns.set(style='white', font_scale=1.8, palette='deep')
        
        # Load the default SPS model.
        t0 = time()
        print('Note: hard-coding zcontinuous!')
        sps = CSPSpecBasis(zcontinuous=1, compute_vega_mags=False)
        print('Initializing the CSPSpecBasis Class took {:.1f} seconds.'.format(time() - t0))

        # Read the parent sample and loop on each object.
        cat = read_parent()
        for obj in cat:
            objprefix = '{0:05}'.format(obj['ISEDFIT_ID'])

            # Grab the emcee / prospector outputs.
            #h5file = os.path.join( datadir(), '{}_{}_mcmc.h5'.format(args.prefix, objprefix) )
            h5file = os.path.join( datadir(), 'test_{}_mcmc.h5'.format(objprefix) )
            print('Reading {}'.format(h5file))

            results, min_results, model = read_results.results_from(h5file, model_file=None)
            
            # Reinitialize the model for this object since it's not written to disk(??).
            model = load_model(results['obs']['zred'])

            # Figure 1: Visualize a random sampling of the MCMC chains.
            chains = np.random.choice(results['run_params']['nwalkers'], size=10, replace=False)
            
            qafile = os.path.join(datadir(), '{}_{}_traces.png'.format(args.prefix, objprefix) )
            print('Generating {}'.format(qafile))
            fig = param_evol(results, figsize=(20, 10), chains=chains)
            #fig.title('Minimization Chains')
            fig.savefig(qafile)

            # Figure 2: Generate a corner/triangle plot of the free parameters.
            params = model.free_params
            nparams = len(params)

            qafile = os.path.join(datadir(), '{}_{}_corner.png'.format(args.prefix, objprefix))
            print('Generating {}'.format(qafile))
            #fig.title('Corners')
            fig = subtriangle(results, start=0, thin=5, truths=None,
                              fig=plt.subplots(nparams, nparams, figsize=(27, 27))[0])
            fig.savefig(qafile)

            # Figure 3: Generate the best-fitting SED.

            # Show the last iteration of a randomly selected walker.
            nwalkers, niter = results['run_params']['nwalkers'], results['run_params']['niter']
            if False:
                theta = results['chain'][nwalkers // 2, niter-1] # initial parameters
            else:
                #print('Plotting based on Powell!!!')
                #theta = min_results.x # initial parameters
                theta = results['model'].initial_theta
            
            mspec, mphot, mextra = model.mean_model(theta, results['obs'], sps=sps)
            print(mextra)

            # Use the filters to set the wavelength and flux limits...
            wspec = sps.csp.wavelengths * (1 + results['obs']['zred']) # spectral wavelengths
            wphot = np.array([f.wave_effective for f in results['obs']['filters']])
            wphot_width = np.array([f.effective_width for f in results['obs']['filters']])

            xmin, xmax = 1000.0, 6E4 
            #xmin, xmax = wphot.min()*0.6, wphot.max()/0.8
            temp = np.interp(np.linspace(xmin, xmax, 10000), wspec, mspec)
            ymin, ymax = -0.1, temp.max()/0.6
            #ymin, ymax = temp.min()*0.8, temp.max()/0.6

            qafile = os.path.join(datadir(), '{}_{}_sed.png'.format(args.prefix, objprefix))
            print('Generating {}'.format(qafile))
            #fig.title('SED Best Fit')
            fig, ax = plt.subplots(figsize=(12, 8))

            # Plot the filter curves...
            if False:
                for ff in range(len(wphot)):
                    f = results['obs']['filters'][ff]
                    w, t = f.wavelength.copy(), f.transmission.copy()
                    while t.max() > 1:
                        t /= 10.0
                        t = 0.1*(ymax-ymin)*t + ymin
                        ax.loglog(w, t, lw=3, color='gray', alpha=0.7)

            wfactor = 1E-4
            factor = 10**(0.4*16.4) # maggies --> mJy
            #factor = 10**(0.4*23.9) # maggies --> microJy
                    
            ax.plot(wfactor * wspec, factor * mspec /  mextra, lw=0.7, alpha=0.7, label='Model spectrum') # color='navy', 
            #ax.loglog(wspec, mspec / mextra, lw=0.7, color='navy', alpha=0.7, label='Model spectrum')
            ax.errorbar(wfactor * wphot, factor * mphot / mextra, marker='s', ls='', lw=3, markersize=20,
                        markerfacecolor='none', markeredgewidth=3, # markeredgecolor='blue', 
                        alpha=0.8, label='Model photometry')
            ax.errorbar(wfactor * wphot, factor * results['obs']['maggies'], yerr=factor * results['obs']['maggies_unc'],
                        marker='o', ls='', lw=3, markersize=10, #markerfacecolor='none',
                        markeredgewidth=3, alpha=0.8, label='Observed photometry')
                        #ecolor='red', markeredgecolor='red', 
                
            ax.set_xlabel(r'Observed-frame Wavelength (${}$m)'.format('\mu'))
            ax.set_ylabel('Flux Density (mJy)')
            ax.set_xlim([wfactor * xmin, wfactor * xmax])
            ax.set_ylim([ymin, factor * ymax])
            ax.legend(loc='upper right', fontsize=20)
            fig.savefig(qafile)

            pdb.set_trace()

if __name__ == "__main__":
    main()
