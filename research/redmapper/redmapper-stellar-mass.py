#!/usr/bin/env python

"""Wrapper on prospector to derive stellar masses for a subset of the redmapper
sample.

python ~/repos/git/siena-astrophysics/research/redmapper/redmapper-stellar-mass.py --verbose --seed 333 --nthreads 24 --dofit > redmapper-100.log 2>& 1 &

"""
from __future__ import division, print_function

import os, sys, pdb
import argparse
import multiprocessing
from time import time, asctime

import numpy as np
import fitsio

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

def read_parent(prefix, first=None, last=None):
    """Read the parent (pilot) catalog."""
    fitsfile = os.path.join(datadir(), '{}_sample.fits'.format(prefix))
    print('Reading {}'.format(fitsfile))
    cat = fitsio.read(fitsfile, ext=1, upper=True)

    if first is not None:
        i0 = first
    else:
        i0 = 0
        
    if last is not None:
        i1 = last + 1
    else:
        i1 = len(cat)

    if i0 == i1:
        indx = np.atleast_1d(i0)
    else:
        indx = np.arange(i0, i1)
    
    return cat, indx

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
    mask = cat['IVARMAGGIES'] > 0
    with np.errstate(invalid='ignore', divide='ignore'):
        obs['maggies_unc'] = np.squeeze(1.0 / np.sqrt(cat['IVARMAGGIES'])) #[:3, :])
    obs['phot_mask'] = mask

    # Input spectroscopy (none for this dataset)
    obs['wavelength'] = None 
    obs['spectrum'] = None
    obs['unc'] = None
    obs['mask'] = None

    # Use initial values based on the iSEDfit results.
    obs['zred'] = cat['Z'] # redshift
    obs['mass'] = 10**cat['MSTAR'] # stellar mass
    obs['logzsol'] = np.log10(cat['ZMETAL'] / 0.019) # stellar metallicity
    obs['tage'] = cat['AGE']  # age
    obs['tau'] = cat['TAU']   # tau (for a delayed SFH)
    obs['dust2'] = 0.1

    # Additional informational keys.
    obs['isedfit_id'] = cat['ISEDFIT_ID']
    
    return obs

def _doleast_squares(lsargs):
    """Fille function for the multiprocessing."""
    return doleast_squares(*lsargs)

def doleast_squares(chivecfn, pinit, chi2args):
    from scipy.optimize import least_squares
    res = least_squares(chivecfn, pinit, method='lm',
                        x_scale='jac', args=chi2args)#, verbose=1)
    return res

def load_model(zred=0.1, seed=None):
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
    from prospector_utilities import logmass2mass

    model_params = []

    ##################################################
    # Fixed priors

    # Galaxy redshift
    model_params.append({
        'name': 'zred',
        'N': 1,
        'isfree': False,
        'init': zred,
        'units': '',
        'prior': None,       
        })

    model_params.append({ # current mass in stars, not integral of SFH
        'name': 'mass_units',
        'N': 1,
        'isfree': False,
        'init': 'mstar', # 'mformed'
        'prior': None,       
        })

    # SFH parameterization (delayed-tau)
    model_params.append({
        'name': 'sfh',
        'N': 1,
        'isfree': False,
        'init':   4, # 4 = delayed tau model
        'units': 'type',
        'prior': None,       
        })

    # IMF (Chabrier)
    model_params.append({
        'name': 'imf_type',
        'N': 1,
        'isfree': False,
        'init':   1, # 1 - Chabrier
        'units': '',
        'prior': None,       
        })

    # Do not include dust emission
    model_params.append({
        'name': 'add_dust_emission',
        'N': 1,
        'isfree': False,
        'init':   False, # do not include dust emission
        'units': 'index',
        'prior': None,       
        })

    ##################################################
    # Free priors / parameters

    # Priors on stellar mass and stellar metallicity
    logmass_prior = priors.TopHat(mini=9.0, maxi=13.0)#, seed=seed)
    logmass_init = np.diff(logmass_prior.range)/2.0 + logmass_prior.range[0] # logmass_prior.sample()
    model_params.append({
        'name': 'logmass',
        'N': 1,
        'isfree': True,
        'init': logmass_init, # mass, 
        'init_disp': 0.5,     # dex
        'units': r'$M_{\odot}$',
        'prior': logmass_prior,
        })

    model_params.append({
        'name': 'mass',
        'N': 1,
        'isfree': False,
        'init': 10**logmass_init,
        'units': '',
        'prior': None,
        'depends_on': logmass2mass,
        })

    logzsol_prior = priors.TopHat(mini=np.log10(0.004/0.019), maxi=np.log10(0.04/0.019))#, seed=seed)
    logzsol_init = np.diff(logzsol_prior.range)/2.0 + logzsol_prior.range[0] # logzsol_prior.sample(), # logzsol,
    model_params.append({
        'name': 'logzsol',
        'N': 1,
        'isfree': True,
        'init': logzsol_init,
        'init_disp': 0.3, # logzsol_prior.range[1] * 0.1,
        'units': r'$\log_{10}\, (Z/Z_\odot)$',
        'prior': logzsol_prior, # roughly (0.2-2)*Z_sun
        })

    # Prior(s) on dust content
    dust2_prior = priors.TopHat(mini=0.0, maxi=3.0)#, seed=seed)
    dust2_init = np.diff(dust2_prior.range)/2.0 + dust2_prior.range[0] # dust2_prior.sample(), # dust2,
    model_params.append({
        'name': 'dust2',
        'N': 1,
        'isfree': True,
        'init': dust2_init,
        'init_disp': 0.5, # dust2_prior.range[1] * 0.1,
        'units': '', # optical depth
        'prior': dust2_prior,
        })
    
    # Priors on tau and age
    #tau_prior = priors.TopHat(mini=0.1, maxi=10.0)#, seed=seed)
    tau_prior = priors.LogUniform(mini=0.1, maxi=10.0)#, seed=seed)
    tau_init = np.diff(tau_prior.range)/2.0 + tau_prior.range[0] # tau_prior.sample(), # tau,
    model_params.append({
        'name': 'tau',
        'N': 1,
        'isfree': True,
        'init': tau_init,
        'init_disp': 1.0, # tau_prior.range[1] * 0.1,
        'units': 'Gyr',
        'prior': tau_prior,
        })

    tage_prior = priors.TopHat(mini=0.5, maxi=15)#, seed=seed)
    tage_init = np.diff(tage_prior.range) / 2.0 + tage_prior.range[0] # tage_prior.sample(), # tage,
    model_params.append( {
        'name': 'tage',
        'N': 1,
        'isfree': True,
        'init': tage_init,
        'init_disp': 2.0, # tage_prior.range[1] * 0.1,
        'units': 'Gyr',
        'prior': tage_prior,
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
    from prospect.likelihood import lnlike_spec, lnlike_phot, chi_spec, chi_phot

    # Calculate the prior probability and exit if not within the prior
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
        'cal': model._speccal, # object calibration spectrum
        }

    # Calculate likelihoods--
    t2 = time()
    lnp_spec = lnlike_spec(model_spec, obs=obs, spec_noise=spec_noise, **vectors)
    lnp_phot = lnlike_phot(model_phot, obs=obs, phot_noise=phot_noise, **vectors)
    d2 = time() - t2
    if False:
        from prospect.likelihood import write_log
        write_log(theta, lnp_prior, lnp_spec, lnp_phot, d1, d2)

    return lnp_prior + lnp_phot + lnp_spec

def chisqfn(theta, model, obs, sps, verbose):
    """Return the negative of lnprobfn for minimization."""
    return -lnprobfn(theta=theta, model=model, obs=obs, sps=sps,
                     verbose=verbose)

def chivecfn(theta, model, obs, sps):
    """Return the residuals instead of a posterior probability or negative chisq,
    for use with least-squares optimization methods.

    """
    resid = lnprobfn(theta=theta, model=model, obs=obs, sps=sps, residuals=True)
    return resid

## MPI pool.  This must be done *after* lnprob and chi2 are defined since slaves
## will only see up to sys.exit()
#try:
#    from emcee.utils import MPIPool
#    pool = MPIPool(debug=False, loadbalance=True)
#    if not pool.is_master():
#        # Wait for instructions from the master process.
#        pool.wait()
#        sys.exit(0)
#except(ImportError, ValueError):
#    pool = None
#    print('Not using MPI.')
#
#def halt(message):
#    """Exit, closing pool safely.
#    """
#    print(message)
#    try:
#        pool.close()
#    except:
#        pass
#    sys.exit(0)

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--prefix', type=str, default='redmapper_sdssphot', help='String to prepend to I/O files.')
    parser.add_argument('--nthreads', type=int, default=16, help='Number of cores to use concurrently.')
    parser.add_argument('--first', type=int, default=None, help='Index of first object to fit.')
    parser.add_argument('--last', type=int, default=None, help='Index of last object to fit.')
    parser.add_argument('--seed', type=int, default=1, help='Random number seed.')
    parser.add_argument('--build-sample', action='store_true', help='Build the sample.')
    parser.add_argument('--dofit', action='store_true', help='Run prospector!')
    parser.add_argument('--refit', action='store_true', help='Refit even if the prospector output files exist.')
    parser.add_argument('--qaplots', action='store_true', help='Make some neat plots.')
    parser.add_argument('--verbose', action='store_true', help='Be loquacious.')
    args = parser.parse_args()

    # Specify the run parameters and initialize the SPS object.
    run_params = {
        'prefix':  args.prefix,
        'verbose': args.verbose,
        'seed':    args.seed,
        # initial optimization choices (nmin is only for L-M optimization)
        'do_levenburg': True,
        'do_powell': False,
        'do_nelder_mead': False,
        'nmin': 10, 
        #'nmin': np.max( (10, args.nthreads) ),
        # emcee fitting parameters
        'nwalkers': 128,
        'nburn': [32, 32, 64], 
        'niter': 256, # 512,
        'interval': 0.1, # save 10% of the chains at a time
        # Nestle fitting parameters
        'nestle_method': 'single',
        'nestle_npoints': 200,
        'nestle_maxcall': int(1e6),
        # Multiprocessing
        'nthreads': args.nthreads,
        # SPS initialization parameters
        'compute_vega_mags': False,
        'vactoair_flag': False, # use wavelengths in air
        'zcontinuous': 1,      # interpolate in metallicity
        }
    print(run_params)
    rand = np.random.RandomState(args.seed)

    if not args.build_sample:
        t0 = time()
        print('Initializing CSPSpecBasis...')
        sps = CSPSpecBasis(zcontinuous=run_params['zcontinuous'],
                           compute_vega_mags=run_params['compute_vega_mags'],
                           vactoair_flag=run_params['vactoair_flag'])
        print('...took {:.1f} seconds.'.format(time() - t0))

    if args.build_sample:
        # Read the parent redmapper catalog, choose a subset of objects and
        # write out.
        cat = read_redmapper()

        # Choose objects with masses from iSEDfit, Kravtsov, and pymorph, but
        # for now just pick three random galaxies.
        #these = np.arange(2) # [300, 301, 302]
        these = np.arange(100) # [300, 301, 302]
        print('Selecting {} galaxies.'.format(len(these)))
        out = cat[these]
        #out = cat[:200]

        outfile = os.path.join(datadir(), '{}_sample.fits'.format(run_params['prefix']))
        print('Writing {}'.format(outfile))
        fitsio.write(outfile, out, clobber=True)

    if args.dofit:
        import h5py
        import emcee
        from prospect import fitting
        from prospect.io import write_results

        # Read the parent sample and loop on each object.
        cat, indx = read_parent(prefix=run_params['prefix'], first=args.first, last=args.last)

        for ii, obj in enumerate(cat[indx]):
            objprefix = '{0:05}'.format(obj['ISEDFIT_ID'])
            print('Working on index {} (fit {}/{} among {} galaxies) with prefix {}.'.format(
                indx[ii], ii+1, len(indx), len(cat), objprefix))

            # Check for the HDF5 output file / fitting results -- 
            outroot = os.path.join( datadir(), '{}_{}'.format(run_params['prefix'], objprefix) )
            hfilename = os.path.join( datadir(), '{}_{}_mcmc.h5'.format(
                run_params['prefix'], objprefix) )
            if os.path.isfile(hfilename):
                if args.refit:
                    os.remove(hfilename)
                else:
                    print('Prospector fitting results {} exist; skipping.'.format(hfilename))
                    continue
            
            # Grab the photometry for this object and then initialize the priors
            # and the SED model.
            obs = getobs(obj)
            model = load_model(zred=obs['zred'], seed=args.seed)

            #pdb.set_trace()
            #sps.csp.params.ssp_params

            # Get close to the right answer doing a simple minimization.
            if run_params['verbose']:
                print('Free parameters: {}'.format(model.free_params))
                print('Initial parameter values: {}'.format(model.initial_theta))
            initial_theta = model.rectify_theta(model.initial_theta) # make zeros tiny numbers

            if bool(run_params.get('do_powell', True)):
                tstart = time()
                # optimization options
                powell_opt = {'ftol': run_params.get('ftol', 0.5e-5),
                              'xtol': 1e-6,
                              'maxfev': run_params.get('maxfev', 5000)}

                chi2args = (model, obs, sps, run_params['verbose']) # extra arguments for chisqfn
                guesses, pinit = fitting.pminimize(chisqfn, initial_theta, args=chi2args, model=model,
                                                   method='powell', opts=powell_opt, pool=pool,
                                                   nthreads=run_params['nthreads'])
                best = np.argmin([p.fun for p in guesses])

                # Hack to recenter values outside the parameter bounds!
                initial_center = fitting.reinitialize(guesses[best].x, model,
                                                      edge_trunc=run_params.get('edge_trunc', 0.1))
                initial_prob = -1 * guesses[best]['fun']
                pdur = time() - tstart
                if run_params['verbose']:
                    print('Powell initialization took {:.1f} seconds.'.format(pdur))
                    print('Best Powell guesses: {}'.format(initial_center))
                    print('Initial probability: {}'.format(initial_prob))
        
            elif bool(run_params.get('do_nelder_mead', True)):
                from scipy.optimize import minimize
                tstart = time()
                chi2args = (model, obs, sps, run_params['verbose']) # extra arguments for chisqfn
                guesses = minimize(chisqfn, initial_theta, args=chi2args, method='nelder-mead')
                pdur = time() - tstart

                # Hack to recenter values outside the parameter bounds!
                initial_center = fitting.reinitialize(guesses.x, model,
                                                      edge_trunc=run_params.get('edge_trunc', 0.1))
                initial_prob = -1 * guesses['fun']
                
                if run_params['verbose']:
                    print('Nelder-Mead initialization took {:.1f} seconds.'.format(pdur))
                    print('Best guesses: {}'.format(initial_center))
                    print('Initial probability: {}'.format(initial_prob))
                    
            elif bool(run_params.get('do_levenburg', True)):
                tstart = time()
                nmin = run_params['nmin']
                
                chi2args = (model, obs, sps)
                pinitial = fitting.minimizer_ball(initial_theta, nmin, model, seed=run_params['seed'])
                
                lsargs = list()
                for pinit in pinitial:
                    lsargs.append((chivecfn, pinit, chi2args))

                if False:
                #if run_params['nthreads'] > 1:
                    p = multiprocessing.Pool(run_params['nthreads'])
                    guesses = p.map(_doleast_squares, lsargs)
                    p.close()
                else:
                    guesses = list()
                    for lsargs1 in lsargs:
                        guesses.append(_doleast_squares(lsargs1))
                
                chisq = [np.sum(r.fun**2) for r in guesses]
                best = np.argmin(chisq)
                initial_prob = -np.log(chisq[best] / 2)
                initial_center = guesses[best].x
                
                pdur = time() - tstart
                if run_params['verbose']:
                    print('Levenburg-Marquardt initialization took {:.1f} seconds.'.format(pdur))
                    print('Best guesses: {}'.format(initial_center))
                    print('Initial probability: {}'.format(initial_prob))

                if False:
                    from prospector_utilities import bestfit_sed
                    fig = bestfit_sed(obs, theta=initial_center, sps=sps, model=model)
                    fig.savefig('test.png')
                    
            else:
                if run_params['verbose']:
                    print('Skipping initial minimization.')
                guesses = None
                pdur = 0.0
                initial_center = initial_theta.copy()
                initial_prob = None

            # Write some basic info to the HDF5 file--
            hfile = h5py.File(hfilename, 'a')
            write_results.write_h5_header(hfile, run_params, model)
            write_results.write_obs_to_h5(hfile, obs)
            
            if run_params['verbose']:
                print('Started emcee sampling on {}'.format(asctime()))
            tstart = time()
            out = fitting.run_emcee_sampler(lnprobfn, initial_center, model,
                                            verbose=run_params['verbose'],
                                            nthreads=run_params['nthreads'],
                                            nwalkers=run_params['nwalkers'],
                                            nburn=run_params['nburn'],
                                            niter=run_params['niter'], 
                                            prob0=initial_prob, hdf5=hfile,
                                            postargs=(model, obs, sps),
                                            pool=None)
            esampler, burn_p0, burn_prob0 = out
            del out
            
            edur = time() - tstart
            if run_params['verbose']:
                print('Finished emcee sampling in {:.2f} minutes.'.format(edur / 60.0))

            # Update the HDF5 file with the results.
            write_results.write_pickles(run_params, model, obs, esampler, guesses,
                                        outroot=outroot, toptimize=pdur, tsample=edur,
                                        sampling_initial_center=initial_center,
                                        post_burnin_center=burn_p0,
                                        post_burnin_prob=burn_prob0)
            write_results.write_hdf5(hfilename, run_params, model, obs, esampler, 
                                     guesses, toptimize=pdur, tsample=edur,
                                     sampling_initial_center=initial_center,
                                     post_burnin_center=burn_p0,
                                     post_burnin_prob=burn_prob0)
            hfile.close()
            
    if args.qaplots:        
        import h5py
        from prospect.io import read_results
        from prospector_utilities import param_evol, subtriangle, bestfit_sed

        # Read the parent sample and loop on each object.
        cat, indx = read_parent(prefix=run_params['prefix'], first=args.first, last=args.last)
        for ii, obj in enumerate(cat[indx]):
            objprefix = '{0:05}'.format(obj['ISEDFIT_ID'])
            
            print('Working on index {} (fit {}/{} among {} galaxies) with prefix {}.'.format(
                indx[ii], ii+1, len(indx), len(cat), objprefix))
    
            # Grab the emcee / prospector outputs.
            h5file = os.path.join( datadir(), '{}_{}_mcmc.h5'.format(run_params['prefix'], objprefix) )
            if not os.path.isfile(h5file):
                print('HDF5 file {} not found; skipping.'.format(h5file))
                continue
            print('Reading {}'.format(h5file))

            results, guesses, model = read_results.results_from(h5file)
            nwalkers, niter, nparams = results['chain'][:, :, :].shape

            # --------------------------------------------------
            # Figure: Generate the best-fitting SED.
            qafile = os.path.join(datadir(), '{}_{}_sed.png'.format(args.prefix, objprefix))
            print('Writing {}'.format(qafile))

            fig = bestfit_sed(results['obs'], chain=results['chain'], lnprobability=results['lnprobability'],
                              sps=sps, model=model, seed=results['run_params']['seed'])
            fig.savefig(qafile)

            # --------------------------------------------------
            # Figure: Visualize a random sampling of the MCMC chains.
            qafile = os.path.join(datadir(), '{}_{}_chains.png'.format(args.prefix, objprefix) )
            print('Writing {}'.format(qafile))

            thesechains = rand.choice(nwalkers, size=int(0.3*nwalkers), replace=False)
            fig = param_evol(results, chains=thesechains)
            fig.savefig(qafile)

            # --------------------------------------------------
            # Figure: Generate a corner/triangle plot of the free parameters.
            qafile = os.path.join(datadir(), '{}_{}_corner.png'.format(args.prefix, objprefix))
            print('Writing {}'.format(qafile))

            fig = subtriangle(results, thin=2)
            fig.savefig(qafile)

if __name__ == "__main__":
    main()
