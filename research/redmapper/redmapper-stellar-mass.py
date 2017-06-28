#!/usr/bin/env python

"""Wrapper on prospector to derive stellar masses for a subset of the redmapper
sample.

"""
from __future__ import division, print_function

import os
import sys
import pdb
import argparse
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
      zred (float): input (fixed) galaxy redshift.
      mass (float): initial stellar mass (Msun)
      logzsol (float): initial metallicity
      tage (float): initial age (Gyr)
      tau (float): initial SFH timescale (Gyr)
      dust2 (float): initial diffuse dust attenuation (dimensionless optical depth)

    Returns:
      sed (prospect.models.sedmodel.SedModel): SED priors and other stuff.

    Notes:
      FSPS parameters are documented here:
        http://dan.iel.fm/python-fsps/current/stellarpop_api/#api-reference

      Initialization parameters:
        * compute_vega_mags (must be set at initialization)
        * vactoair_flag (must be set at initialization)
        * zcontinuous (must be set at initialization)
    
      Metallicity parameters:
        * zmet (default 1, ignored if zcontinuous>0)
        * logzsol (default 0.0, used if zcontinuous>0)
        * pmetals (default 2.0, only used if zcontinuous=2)

      Dust parameters:
        * add_agb_dust_model (default True)
        * add_dust_emission (default True)
        * cloudy_dust (default False)
        * agb_dust (default 1.0)
        * dust_type (default 0=power law)
        * dust_index, dust1_index
        * dust_tesc
        * dust1 (default 0.0) - extra optical depth towards young stars at 5500A
        * dust2 (default 0.0) - diffuse dust optical depth towards all stars at 5500A
        * dust_clumps, frac_nodust, frac_obrun
        * mwr, uvb, wgp1, wgp2, wgp3, 
        * duste_gamma, duste_umin, duste_qpah

      Star formation history parameters:
        * sfh (default 0=SSP, 1=tau, 4=delayed, 5=truncated delayed tau)
        * tau (default 1)
        * const, sf_start, sf_trunc
        * tage (default 0.0)
        * fburst, tburst, sf_slope
    
      Miscellaneous parameters:
        * add_igm_absorption (default False)
        * igm_factor (default 1.0)
        * smooth_velocity (default True)
        * sigma_smooth, min_wave_smooth, max_wave_smooth
        * redshift_colors (default False, do not use)
        * compute_light_ages (default False, do not use)
       
      Stellar population parameters:
        * add_stellar_remnants (default True)
        * tpagb_norm_type (default 2)
        * dell (default 0.0, do not use)
        * delt (default 0.0, do not use)
        * redgb (default 1.0)
        * fcstar (default 1.0)
        * sbss (default 0.0)
        * fbhb (default 0.0)
        * pagb (default 1.0)
        * imf_type (default 2=Kroupa01)
        * imf1, imf2, imf3, vdmc, mdave, masscut
        * evtype (default 1)
        * tpagb_norm_type

      Emission lines:
        * add_neb_emission (default False)
        * add_neb_continuum (default False)
        * gas_logz (default 0.0)
        * gas_logu (default -2)

      Galaxy properties:
        * zred (default 0.0)

      AGN properties:
        * fagn (default 0.0)
        * agn_tau (default 10)

      Calibration parameters:
        * phot_jitter

    """
    from prospect.models import priors, sedmodel
    
    model_params = []

    # (Fixed) prior on galaxy redshift.
    model_params.append({
        'name': 'zred',
        'N': 1,
        'isfree': False,
        'init': zred,
        'units': '',
        })
    
    # Priors on stellar mass and stellar metallicity.
    model_params.append({
        'name': 'mass_units',
        'N': 1,
        'isfree': False,
        'init': 'mformed'
        })

    model_params.append({
        'name': 'mass',
        'N': 1,
        'isfree': True,
        'init':      mass, 
        'init_disp': 5e11, 
        'units': r'$M_{\odot}$',
        'prior_function': priors.LogUniform(mini=1e8, maxi=1e13),
        })

    model_params.append({
        'name': 'logzsol',
        'N': 1,
        'isfree': False,
        'init': logzsol,
        'init_disp': 0.3, # dex
        'units': r'$\log_{10}\, (Z/Z_\odot)$',
        'prior_function': priors.TopHat(mini=-1.5, maxi=0.19),
        })

    # Priors on dust
    model_params.append({
        'name': 'dust2',
        'N': 1,
        'isfree': False,
        'init': dust2,
        'init_disp': 0.3,
        'units': '',
        'prior_function': priors.TopHat(mini=0.0, maxi=2.0),
        })
    
    # Prior on the IMF.
    model_params.append({
        'name': 'imf_type',
        'N': 1,
        'isfree': False,
        'init':   1, # 1 - Chabrier
        'units': ''
        })

    # Priors on SFH type (fixed), tau, and age.
    model_params.append({
        'name': 'sfh',
        'N': 1,
        'isfree': False,
        'init':   4, # 4 = delayed tau model
        'units': 'type'
        })

    model_params.append({
        'name': 'tau',
        'N': 1,
        'isfree': True,
        'init': tau,
        'init_disp': 1.0,
        'units': 'Gyr',
        'prior_function': priors.LogUniform(mini=0.1, maxi=10.0),
        })

    model_params.append( {
        'name':   'tage',
        'N':       1,
        'isfree':    True,
        'init':      tage,
        'init_disp':  3.0,
        'units':       'Gyr',
        'prior_function': priors.TopHat(mini=0.5, maxi=15),
        })

    model = sedmodel.SedModel(model_params)
    
    return model

def lnprobfn(theta, model, obs, sps, verbose=False, spec_noise=None,
             phot_noise=None, residuals=False):
    """Define the likelihood function.

    Given a parameter vector and a dictionary of observational data and a model
    object, return the ln of the posterior. This requires that an sps object
    (and if using spectra and gaussian processes, a GP object) be instantiated.

    """
    from prospect.likelihood import (lnlike_spec, lnlike_phot, write_log,
                                     chi_spec, chi_phot)

    lnp_prior = model.prior_product(theta)
    if not np.isfinite(lnp_prior):
        return -np.infty
        
    # Generate the mean model--
    t1 = time()
    try:
        model_spec, model_phot, model_extras = model.mean_model(theta, obs, sps=sps)
    except(ValueError):
        return -np.infty
    d1 = time() - t1

    # Return chi vectors for least-squares optimization--
    if residuals:
        chispec = chi_spec(model_spec, obs)
        chiphot = chi_phot(model_phot, obs)
        return np.concatenate([chispec, chiphot])

    # Noise modeling--
    if spec_noise:
        spec_noise.update(**model.params)
    if phot_noise:
        phot_noise.update(**model.params)

    vectors = {
        'spec': model_spec,    # model spectrum
        'phot': model_phot,    # model photometry
        'sed': model._spec,    # object spectrum
        #'unc': obs['unc'],     # object uncertainty spectrum
        'cal': model._speccal, # object calibration spectrum
        #'maggies_unc': obs['maggies_unc'] # object photometric uncertainty
        }

    # Calculate log-likelihoods--
    t2 = time()
    lnp_spec = lnlike_spec(model_spec, obs=obs, spec_noise=spec_noise, **vectors)
    lnp_phot = lnlike_phot(model_phot, obs=obs, phot_noise=phot_noise, **vectors)
    d2 = time() - t2
    if verbose:
        write_log(theta, lnp_prior, lnp_spec, lnp_phot, d1, d2)

    return lnp_prior + lnp_phot + lnp_spec

def chisqfn(theta, model, obs, sps, verbose):
    """Return the negative of lnprobfn for minimization."""
    return -lnprobfn(theta, model, obs, sps, verbose)

def chivecfn(theta, model, obs, sps, verbose):
    """Return the residuals instead of a posterior probability or negative
    chisq, for use with least-squares optimization methods
    """
    return lnprobfn(theta, model, obs, sps, verbose, residuals=True)

# MPI pool.  This must be done *after* lnprob and chi2 are defined since slaves
# will only see up to sys.exit()
try:
    from emcee.utils import MPIPool
    pool = MPIPool(debug=False, loadbalance=True)
    if not pool.is_master():
        # Wait for instructions from the master process.
        pool.wait()
        sys.exit(0)
except(ImportError, ValueError):
    pool = None
    print('Not using MPI.')

def halt(message):
    """Exit, closing pool safely.
    """
    print(message)
    try:
        pool.close()
    except:
        pass
    sys.exit(0)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', type=str, default='test', help='String to prepend to I/O files.')
    parser.add_argument('--build-sample', action='store_true', help='Build the sample.')
    parser.add_argument('--do-fit', action='store_true', help='Run prospector!')
    parser.add_argument('--qaplots', action='store_true', help='Make some neat plots.')
    parser.add_argument('--threads', default=16, help='Number of cores to use concurrently.')
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
        from prospect import fitting
        from prospect.io import write_results

        # Specify the run parameters.
        run_params = {
            'prefix':   args.prefix,
            'verbose':  args.verbose,
            'debug':    False,
            # initial optimization choices
            'do_levenburg': True,
            'do_powell': False,
            'do_nelder_mead': False,
            'nmin': 10, # number of starting conditions to sample from the prior for L-M optimization
            # emcee fitting parameters
            'nwalkers': 128,
            'nburn': [32, 32, 64], 
            'niter': 128,
            'interval': 0.1, # save 10% of the chains at a time
            # Nestle fitting parameters
            'nestle_method': 'single',
            'nestle_npoints': 200,
            'nestle_maxcall': int(1e6),
            # SPS initialization parameters
            'compute_vega_mags': False,
            'vactoair_flag':      True, # use wavelengths in air
            'zcontinuous': 1,           # interpolate in metallicity
            }

        # Initialize the SPS object.
        t0 = time()
        print('Initializing the CSPSpecBasis object...')
        sps = CSPSpecBasis(zcontinuous=run_params['zcontinuous'],
                           compute_vega_mags=run_params['compute_vega_mags']) #,
                           #vactoair_flag=run_params['vactoair_flag'])
        print('...took {:.1f} seconds.'.format(time() - t0))
            
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
            if run_params['verbose']:
                print('Initial parameter values: {}'.format(model.initial_theta))
            initial_theta = model.rectify_theta(model.initial_theta) # make zeros tiny numbers

            if bool(run_params.get('do_powell', True)):
                ts = time()
                # optimization options
                powell_opt = {'ftol': run_params.get('ftol', 0.5e-5),
                              'xtol': 1e-6,
                              'maxfev': run_params.get('maxfev', 5000)}

                chi2args = (sps, run_params['verbose'])
                guesses, pinit = fitting.pminimize(chisqfn, initial_theta, args=chi2args, model=model,
                                                   method='powell', opts=powell_opt, pool=pool,
                                                   nthreads=run_params.get('nthreads', 1))
                best = np.argmin([p.fun for p in guesses])

                # Hack to recenter values outside the parameter bounds!
                initial_center = fitting.reinitialize(guesses[best].x, model,
                                                      edge_trunc=run_params.get('edge_trunc', 0.1))
                initial_prob = -1 * guesses[best]['fun']
                pdur = time() - ts
                if run_params['verbose']:
                    print('Powell initialization took {} seconds.'.format(pdur))
                    print('Best Powell guesses: {}'.format(initial_center))
                    print('Initial probability: {}'.format(initial_prob))
        
            elif bool(run_params.get('do_nelder_mead', True)):
                from scipy.optimize import minimize
                ts = time()
                chi2args = (model, obs, sps, run_params['verbose']) # extra arguments for chisqfn
                guesses = minimize(chisqfn, initial_theta, args=chi2args, method='nelder-mead')
                pdur = time() - ts

                # Hack to recenter values outside the parameter bounds!
                initial_center = fitting.reinitialize(guesses.x, model,
                                                      edge_trunc=run_params.get('edge_trunc', 0.1))
                initial_prob = -1 * guesses['fun']
                
                if run_params['verbose']:
                    print('Nelder-Mead initialization took {} seconds.'.format(pdur))
                    print('Best guesses: {}'.format(initial_center))
                    print('Initial probability: {}'.format(initial_prob))
                    
            elif bool(run_params.get('do_levenburg', True)):
                from scipy.optimize import least_squares
                ts = time()
                nmin = run_params.get('nmin', 10)
                
                chi2args = (model, obs, sps, run_params['verbose']) # extra arguments for chisqfn
                pinitial = fitting.minimizer_ball(model.initial_theta.copy(), nmin, model)
                guesses = []
                for i, pinit in enumerate(pinitial):
                    res = least_squares(chivecfn, pinit, method='lm', x_scale='jac', args=chi2args)
                    guesses.append(res)
        
                chisq = [np.sum(r.fun**2) for r in guesses]
                best = np.argmin(chisq)

                # Hack to recenter values outside the parameter bounds!
                initial_center = fitting.reinitialize(guesses[best].x, model,
                                                      edge_trunc=run_params.get('edge_trunc', 0.1))
                initial_prob = None
                pdur = time() - ts
                if run_params['verbose']:
                    print('Levenburg-Marquardt initialization took {} seconds.'.format(pdur))
                    print('Best guesses: {}'.format(initial_center))
                    print('Initial probability: {}'.format(initial_prob))
        
            else:
                if run_params['verbose']:
                    print('Skipping initial minimization.')
                guesses = None
                pdur = 0.0
                initial_center = initial_theta.copy()
                initial_prob = None

            pdb.set_trace()

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
