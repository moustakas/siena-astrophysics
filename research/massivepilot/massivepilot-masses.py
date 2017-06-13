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

##########################################333
# MOVE TO ITS OWN FUNCTION!

import sys, os
import numpy as np

try:
    from sedpy.observate import load_filters
except:
    pass

"""Convenience functions for reading and reconstructing results from a fitting
run, including reconstruction of the model for making posterior samples
"""

#__all__ = ["subtriangle", "param_evol"]

def model_comp(theta, model, obs, sps, photflag=0, gp=None):
    """Generate and return various components of the total model for a given
    set of parameters.
    """
    logarithmic = obs.get('logify_spectrum')
    obs, _, _ = obsdict(obs, photflag=photflag)
    mask = obs['mask']
    mu = model.mean_model(theta, obs, sps=sps)[photflag][mask]
    sed = model.sed(theta, obs, sps=sps)[photflag][mask]
    wave = obs['wavelength'][mask]

    if photflag == 0:
        cal = model.spec_calibration(theta, obs)
        if type(cal) is not float:
            cal = cal[mask]
        try:
            s, a, l = model.spec_gp_params()
            gp.kernel[:] = np.log(np.array([s[0], a[0]**2, l[0]**2]))
            spec = obs['spectrum'][mask]
            if logarithmic:
                gp.compute(wave, obs['unc'][mask])
                delta = gp.predict(spec - mu, wave)
            else:
                gp.compute(wave, obs['unc'][mask], flux=mu)
                delta = gp.predict(spec - mu)
            if len(delta) == 2:
                delta = delta[0]
        except(TypeError, AttributeError, KeyError):
            delta = 0
    else:
        mask = np.ones(len(obs['wavelength']), dtype=bool)
        cal = np.ones(len(obs['wavelength']))
        delta = np.zeros(len(obs['wavelength']))

    return sed, cal, delta, mask, wave


def param_evol(sample_results, showpars=None, start=0, figsize=None, chains=None, **plot_kwargs):
    """Plot the evolution of each parameter value with iteration #, for each
    walker in the chain.

    :param sample_results:
        A Prospector results dictionary, usually the output of
        ``results_from('resultfile')``.

    :param showpars: (optional)
        A list of strings of the parameters to show.  Defaults to all
        parameters in the ``"theta_labels"`` key of the ``sample_results``
        dictionary.

    :param start: (optional, default: 0)
        Integer giving the iteration number from which to start plotting.

    :param **plot_kwargs:
        Extra keywords are passed to the
        ``matplotlib.axes._subplots.AxesSubplot.plot()`` method.

    :returns tracefig:
        A multipaneled Figure object that shows the evolution of walker
        positions in the parameters given by ``showpars``, as well as
        ln(posterior probability)
    """
    import matplotlib.pyplot as pl

    if chains is None:
        chain = sample_results['chain'][:, start:, :]
        lnprob = sample_results['lnprobability'][:, start:]
    else:
        chain = sample_results['chain'][chains, start:, :]
        lnprob = sample_results['lnprobability'][:, start:]
    nwalk = chain.shape[0]
    try:
        parnames = np.array(sample_results['theta_labels'])
    except(KeyError):
        print('FIX THIS!!!!!!!!!!!!!!!!!!!!!!!!!!')
        parnames = np.array(['Bob1', 'Bob2', 'Bob3'])
        #parnames = np.array(sample_results['model'].theta_labels())

    # logify mass
    if 'mass' in parnames:
        midx = [l=='mass' for l in parnames]
        chain[:,:,midx] = np.log10(chain[:,:,midx])
        parnames[midx] = 'logmass'

    # Restrict to desired parameters
    if showpars is not None:
        ind_show = np.array([p in showpars for p in parnames], dtype=bool)
        parnames = parnames[ind_show]
        chain = chain[:, :, ind_show]

    # Set up plot windows
    ndim = len(parnames) + 1
    nx = int(np.floor(np.sqrt(ndim)))
    ny = int(np.ceil(ndim * 1.0 / nx))
    sz = np.array([nx, ny])
    factor = 3.0           # size of one side of one panel
    lbdim = 0.2 * factor   # size of left/bottom margin
    trdim = 0.2 * factor   # size of top/right margin
    whspace = 0.05 * factor         # w/hspace size
    plotdim = factor * sz + factor * (sz - 1) * whspace
    dim = lbdim + plotdim + trdim

    if figsize is None:
        fig, axes = pl.subplots(nx, ny, figsize=(dim[1], dim[0]))
    else:
        fig, axes = pl.subplots(nx, ny, figsize=figsize)
    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(left=lb[1], bottom=lb[0], right=tr[1], top=tr[0],
                        wspace=whspace, hspace=whspace)

    # Sequentially plot the chains in each parameter
    for i in range(ndim - 1):
        ax = axes.flatten()[i]
        for j in range(nwalk):
            ax.plot(chain[j, :, i], **plot_kwargs)
        ax.set_title(parnames[i], y=1.02)
    # Plot lnprob
    ax = axes.flatten()[-1]
    for j in range(nwalk):
        ax.plot(lnprob[j, :])
    ax.set_title('lnP', y=1.02)
    pl.tight_layout()
    return fig


def subtriangle(sample_results, outname=None, showpars=None,
                start=0, thin=1, truths=None, trim_outliers=None,
                extents=None, **kwargs):
    """Make a triangle plot of the (thinned, latter) samples of the posterior
    parameter space.  Optionally make the plot only for a supplied subset of
    the parameters.

    :param start:
        The iteration number to start with when drawing samples to plot.

    :param thin:
        The thinning of each chain to perform when drawing samples to plot.

    :param showpars:
        List of string names of parameters to include in the corner plot.

    :param truths:
        List of truth values for the chosen parameters
    """
    try:
        import triangle
    except(ImportError):
        import corner as triangle

    # pull out the parameter names and flatten the thinned chains
    try:
        parnames = np.array(sample_results['theta_labels'])
    except(KeyError):
        #parnames = np.array(sample_results['model'].theta_labels())
        print('FIX THIS!!!!!!!!!!!!!!!!!!!!!!!!!!')
        parnames = np.array(['Bob1', 'Bob2', 'Bob3'])

    flatchain = sample_results['chain'][:, start::thin, :]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])

    # logify mass
    if 'mass' in parnames:
        midx = [l=='mass' for l in parnames]
        flatchain[:,midx] = np.log10(flatchain[:,midx])
        parnames[midx] = 'logmass'

    # restrict to parameters you want to show
    if showpars is not None:
        ind_show = np.array([p in showpars for p in parnames], dtype=bool)
        flatchain = flatchain[:, ind_show]
        #truths = truths[ind_show]
        parnames = parnames[ind_show]
    if trim_outliers is not None:
        trim_outliers = len(parnames) * [trim_outliers]
    try:
        fig = triangle.corner(flatchain, labels=parnames, truths=truths,  verbose=False,
                              quantiles=[0.16, 0.5, 0.84], range=trim_outliers, **kwargs)
    except:
        fig = triangle.corner(flatchain, labels=parnames, truths=truths,  verbose=False,
                              quantiles=[0.16, 0.5, 0.84], range=trim_outliers, **kwargs)

    if outname is not None:
        fig.savefig('{0}.triangle.png'.format(outname))
        #pl.close(fig)
    else:
        return fig


def obsdict(inobs, photflag):
    """Return a dictionary of observational data, generated depending on
    whether you're matching photometry or spectroscopy.
    """
    obs = inobs.copy()
    if photflag == 0:
        outn = 'spectrum'
        marker = None
    elif photflag == 1:
        outn = 'sed'
        marker = 'o'
        obs['wavelength'] = np.array([f.wave_effective for f in obs['filters']])
        obs['spectrum'] = obs['maggies']
        obs['unc'] = obs['maggies_unc']
        obs['mask'] = obs['phot_mask'] > 0

    return obs, outn, marker

# All this because scipy changed the name of one class, which
# shouldn't even be a class.

renametable = {
    'Result': 'OptimizeResult',
    }


def mapname(name):
    if name in renametable:
        return renametable[name]
    return name


def mapped_load_global(self):
    module = mapname(self.readline()[:-1])
    name = mapname(self.readline()[:-1])
    klass = self.find_class(module, name)
    self.append(klass)


def load(file):
    unpickler = pickle.Unpickler(file)
    unpickler.dispatch[pickle.GLOBAL] = mapped_load_global
    return unpickler.load()

##########################################333
##########################################333
##########################################333
##########################################333
##########################################333

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

    # Input redshift.
    obs['zred'] = cat['Z']

    # Additional informational keys.
    obs['isedfit_id'] = cat['ISEDFIT_ID']
    
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
                         'init': 5e11, # changed from 1-5 06/07/17 # C
                         'init_disp': 5e10, # changed from 1-5 06/07/17 # C 
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
                         'init_disp': 1.5, # changed from 0.5-1.5  06/07/17 # C 
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
        these = [3, 4, 5]
        print('Selecting {} galaxies.'.format(len(these)))
        out = cat[these]

        outfile = os.path.join(datadir(), 'massivepilot-sample.fits')
        print('Writing {}'.format(outfile))
        fitsio.write(outfile, out)
        
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
            objprefix = '{0:05}'.format(obj['ISEDFIT_ID'])
            
            # Grab the photometry for this object and then initialize the priors
            # and the SED model.
            obs = getobs(obj)
            model = load_model(obs['zred'])

            # Get close to the right answer doing a simple minimization.
            initial_theta = model.rectify_theta(model.initial_theta) # initial parameters
            
            t0 = time()
            min_results = minimize(chisqfn, initial_theta, (model, obs, sps), method=min_method)
            pdur = time() - t0

            print('What is edge_trunc!')
            initial_center = fitting.reinitialize(
                min_results.x, model, edge_trunc=run_params.get('edge_trunc', 0.1)
                ) #uncommented edge_trunc 06/07/17 # C
            initial_prob = -1 * min_results['fun']
            
            print('Minimization {} finished in {} seconds'.format(min_method, pdur))
            print('best {0} guess: {1}'.format(min_method, initial_center))
            print('best {0} lnp: {1}'.format(min_method, initial_prob))
            
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
            out = fitting.run_emcee_sampler(lnprobfn, initial_center, model,
                                            threads=args.threads, 
                                            initial_prob=initial_prob,
                                            hdf5=hfile, nwalkers=run_params.get('nwalkers'),
                                            nburn=run_params.get('nburn'),
                                            niter=run_params.get('niter'), 
                                            postargs=(model, obs, sps))
            esampler, burn_p0, burn_prob0 = out
            edur = time() - tstart

            sys.stdout = fout

            print('done emcee in {}s'.format(edur))
            
            # Write out more.
        
            write_results.write_pickles(run_params, model, obs, esampler, min_results,
                                        outroot=run_params['prefix'], toptimize=pdur, tsample=edur,
                                        sampling_initial_center=initial_center,
                                        post_burnin_center=burn_p0,
                                        post_burnin_prob=burn_prob0)
            pdb.set_trace()
            write_results.write_hdf5(hfilename, run_params, model, esampler, min_results,
                                     min_results, toptimize=pdur, tsample=edur,
                                     sampling_initial_center=initial_center,
                                     post_burnin_center=burn_p0,
                                     post_burnin_prob=burn_prob0)

    if args.qaplots:
        from prospect.io import read_results

        # Load the default SPS model.
        t0 = time()
        print('Do we really need to do this here???  Note: hard-coding zcontinuous!')
        sps = CSPSpecBasis(zcontinuous=1, compute_vega_mags=False)
        print('Initializing the CSPSpecBasis Class took {:.1f} seconds.'.format(time() - t0))

        # Read the parent sample and loop on each object.
        cat = read_parent()
        for obj in cat:
            objprefix = '{0:05}'.format(obj['ISEDFIT_ID'])

            # Grab the emcee / prospector outputs.
            h5file = os.path.join( datadir(), '{}_{}_mcmc.h5'.format(args.prefix, objprefix) )
            print('Reading {}'.format(h5file))

            results, powell_results, model = read_results.results_from(h5file)

            # Reinitialize the model for this object since it's not written to disk(??).
            model = load_model(results['obs']['zred'])

            # Figure 1: Visualize a random sampling of the MCMC chains.
            chains = np.random.choice(results['run_params']['nwalkers'], size=10, replace=False)
            
            qafile = os.path.join(datadir(), '{}_{}_traces.png'.format(args.prefix, objprefix) )
            print('Generating {}'.format(qafile))
            fig = param_evol(results, figsize=(20, 10), chains=chains)
            fig.savefig(qafile)

            # Figure 2: Generate a corner/triangle plot of the free parameters.
            params = model.free_params
            nparams = len(params)
            #theta_truth = np.array([results['run_params'][pp] for pp in params])

            qafile = os.path.join(datadir(), '{}_{}_corner.png'.format(args.prefix, objprefix))
            print('Generating {}'.format(qafile))
            fig = subtriangle(results, start=0, thin=5, truths=None,
                              fig=plt.subplots(nparams, nparams, figsize=(27, 27))[0])
            fig.savefig(qafile)

            # Figure 3: Generate the best-fitting SED.

            # Show the last iteration of a randomly selected walker.
            nwalkers, niter = results['run_params']['nwalkers'], results['run_params']['niter']
            theta = results['chain'][nwalkers // 2, niter-1] # initial parameters

            mspec, mphot, mextra = model.mean_model(theta, results['obs'], sps=sps)

            # Use the filters to set the wavelength and flux limits...
            wspec = sps.csp.wavelengths # spectral wavelengths
            wphot = np.array([f.wave_effective for f in results['obs']['filters']])
            wphot_width = np.array([f.effective_width for f in results['obs']['filters']])

            xmin, xmax = wphot.min()*0.8, wphot.max()/0.8
            temp = np.interp(np.linspace(xmin, xmax, 10000), wspec, mspec)
            ymin, ymax = temp.min()*0.8, temp.max()/0.8

            qafile = os.path.join(datadir(), '{}_{}_sed.png'.format(args.prefix, objprefix))
            print('Generating {}'.format(qafile))
            fig, ax = plt.subplots(figsize=(16, 8))

            # Plot the filter curves...
            for ff in range(len(wphot)):
                f = results['obs']['filters'][ff]
                w, t = f.wavelength.copy(), f.transmission.copy()
                while t.max() > 1:
                    t /= 10.0
                    t = 0.1*(ymax-ymin)*t + ymin
                    ax.loglog(w, t, lw=3, color='gray', alpha=0.7)
                    
            ax.loglog(wspec, mspec, lw=0.7, color='navy', alpha=0.7, label='Model spectrum')
            ax.errorbar(wphot, mphot, marker='s', ls='', lw=3, markersize=10, markerfacecolor='none',
                        markeredgecolor='blue', markeredgewidth=3, alpha=0.8, label='Model photometry')
            ax.errorbar(wphot, results['obs']['maggies'], yerr=results['obs']['maggies_unc'],
                        ecolor='red', marker='o', ls='', lw=3, markersize=10, markerfacecolor='none',
                        markeredgecolor='red', markeredgewidth=3,
                        alpha=0.8, label='Observed photometry')
                
            ax.set_xlabel('Wavelength [A]')
            ax.set_ylabel('Flux Density [maggies]')
            ax.set_xlim([xmin, xmax])
            ax.set_ylim([ymin, ymax])
            ax.legend(loc='lower right', fontsize=20)
            fig.savefig(qafile)

            pdb.set_trace()

if __name__ == "__main__":
    main()
