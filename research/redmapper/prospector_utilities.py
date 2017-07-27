"""Convenience functions for reading and reconstructing results from a fitting
run, including reconstruction of the model for making posterior samples.

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='white', font_scale=1.7, palette='deep')
setcolors = sns.color_palette()
    
#sns.set_style({'axes.linewidth' : 3.0, 'lines.linewidth': 3, 'lines.markersize': 30})
#sns.set(style='white', font_scale=1.8, palette='deep')
#sns.set(style='white', font_scale=3, palette='deep')
#sns.set(style='white', font_scale=4, palette='dark')
#sns.set_context("talk", rc={"lines.linewidth": 10})
#sns.set(style='white', font_scale=1.8, palette='deep')

ang2micron = 1e-4 # Angstrom --> micron
maggies2mJy = 10**(0.4*16.4) # maggies --> mJy
#maggies2muJy = 10**(0.4*23.9) # maggies --> microJy

def logmass2mass(logmass=11.0, **extras):
    return 10**logmass

def datadir():
    """Return the top-level data directory."""
    return os.path.join(os.getenv('HOME'), 'stellar-mass')    

def _niceparnames(parnames):
    """Replace parameter names with nice names."""

    old = list(['tau',
           'tage',
           'mass',
           'logmass',
           'logzsol',
           'dust2'])
    new = list([r'$\tau$ (Gyr)',
           'Age (Gyr)',
           r'$M / M_{\odot}$',
           r'$\log_{10}\,(M / M_{\odot})$',
           r'$\log_{10}\, (Z / Z_{\odot})$',
           r'$\tau_{diffuse}$'])

    niceparnames = list(parnames).copy()
    for oo, nn in zip( old, new ):
        this = np.where(np.in1d(parnames, oo))[0]
        if len(this) > 0:
            niceparnames[this[0]] = nn
            
    return np.array(niceparnames)

def _galaxyphot(obs):
    """Get the galaxy photometry and inverse variances (converted to mJy) and filter
    effective wavelengths (converted to microns).

    """
    weff = np.array([f.wave_effective for f in obs['filters']]) * ang2micron
    fwhm = np.array([f.effective_width for f in obs['filters']]) * ang2micron

    galphot = obs['maggies'] * maggies2mJy
    galphoterr = obs['maggies_unc'] * maggies2mJy

    return weff, fwhm, galphot, galphoterr

def _sed(model, theta, obs, sps):
    """Construct the SED for a given set of parameters.  Divide by mextra to account
    for the *current* mass in stars (rather than the integrated stellar mass
    based on the SFH.

    Also convert wavelengths from Angstroms to microns and fluxes from maggies
    to mJy.

    """
    modelwave = sps.csp.wavelengths * (1 + obs['zred']) # [observed-frame wavelengths]
    modelwave *= ang2micron 
    
    modelspec, modelphot, mextra = model.mean_model(theta, obs, sps=sps)
    modelspec = modelspec * maggies2mJy
    modelphot = modelphot * maggies2mJy
    
    return modelwave, modelspec, modelphot

def bestfit_sed(obs, chain=None, lnprobability=None, theta=None, sps=None,
                model=None, seed=None, nrand=100):
    """Plot the (photometric) best-fitting SED.

    Either pass chain and lnprobability (to visualize the emcee fitting results)
    *or* theta (to visualize just a single SED fit).

    """
    rand = np.random.RandomState(seed)

    # Get the galaxy photometry and filter info.
    weff, fwhm, galphot, galphoterr = _galaxyphot(obs)

    # Build the maximum likelihood model fit and also grab a random sampling of
    # the chains with weight equal to the posterior probability.    
    if chain is not None:
        nwalkers, niter, nparams = chain.shape
        ntot = nwalkers * niter

        flatchain = chain.reshape(ntot, nparams)
        lnp = lnprobability.reshape(ntot)

        theta = flatchain[lnp.argmax(), :] # maximum likelihood values

        prob = np.exp(lnp - lnp.max())
        prob /= prob.sum()
        rand_indx = rand.choice(ntot, size=nrand, replace=False, p=prob)
        theta_rand = flatchain[rand_indx, :]
        
    modelwave, modelspec, modelphot = _sed(model=model, theta=theta, obs=obs, sps=sps)

    # Establish the wavelength and flux limits.
    #minwave, maxwave = 0.1, 6.0
    minwave, maxwave = np.min(weff - 5*fwhm), np.max(weff + fwhm)

    inrange = (modelwave > minwave) * (modelwave < maxwave)
    maxflux = np.hstack( (galphot + 3*galphoterr, modelspec[inrange]) ).max() * 1.05
    minflux = -0.05 * maxflux

    fig, ax = plt.subplots(figsize=(12, 8))
    if chain is not None:
        for ii in range(nrand):
            _, r_modelspec, _ = _sed(model=model, theta=theta_rand[ii, :], obs=obs, sps=sps)
            ax.plot(modelwave, r_modelspec, alpha=0.2, color='gray')
    ax.plot(modelwave, modelspec, alpha=1.0, label='Model spectrum')
    
    ax.errorbar(weff, modelphot, marker='s', ls='', lw=3, markersize=20, markerfacecolor='none',
                markeredgewidth=3, alpha=0.8, label='Model photometry')
    ax.errorbar(weff, galphot, yerr=galphoterr, marker='o', ls='', lw=3, markersize=10,
                markeredgewidth=3, alpha=0.8, label='Observed photometry')
                
    ax.set_xlabel(r'Observed-frame Wavelength (${}$m)'.format('\mu'))
    ax.set_ylabel('Flux Density (mJy)')
    ax.set_xlim(minwave, maxwave)
    ax.set_ylim(minflux, maxflux)
    ax.legend(loc='upper right', fontsize=16, frameon=True)
    fig.subplots_adjust(left=0.1, right=0.95, bottom=0.12, top=0.95)

    return fig

def param_evol(sample_results, showpars=None, start=0, figsize=None,
               chains=None, **plot_kwargs):
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
    if chains is None:
        chain = sample_results['chain'][:, start:, :]
        lnprob = sample_results['lnprobability'][:, start:]
    else:
        chain = sample_results['chain'][chains, start:, :]
        lnprob = sample_results['lnprobability'][:, start:]
    nwalk = chain.shape[0]
    
    parnames = sample_results['theta_labels']

    # logify mass
    #if 'mass' in parnames:
    if 'mass' in parnames and 'logmass' not in parnames:
        midx = np.where(np.in1d(parnames, 'mass'))[0]
        if len(midx) > 0:
            chain[:, :, midx[0]] = np.log10(chain[:, :, midx[0]])
            parnames[midx[0]] = 'logmass'
    parnames = np.array(parnames)

    # Restrict to desired parameters
    if showpars is not None:
        ind_show = np.array([p in showpars for p in parnames], dtype=bool)
        parnames = parnames[ind_show]
        chain = chain[:, :, ind_show]

    # Make nice labels.
    niceparnames = _niceparnames(parnames)

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
        fig, axes = plt.subplots(nx, ny, figsize=(dim[1], dim[0]))
    else:
        fig, axes = plt.subplots(nx, ny, figsize=figsize)
    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    fig.subplots_adjust(left=lb[1], bottom=lb[0], right=tr[1], top=tr[0],
                        wspace=whspace, hspace=whspace)

    # Sequentially plot the chains in each parameter
    for i in range(ndim - 1):
        ax = axes.flatten()[i]
        for j in range(nwalk):
            ax.plot(chain[j, :, i], **plot_kwargs)
        ax.set_title(niceparnames[i], y=1.02)
        
    # Plot lnprob
    ax = axes.flatten()[-1]
    for j in range(nwalk):
        ax.plot(lnprob[j, :])
    ax.set_title(r'$\ln\,P$', y=1.02)
    plt.tight_layout()

    return fig

def subtriangle(sample_results, outname=None, showpars=None, start=0, thin=1,
                truths=None, trim_outliers=None, extents=None, **kwargs):
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
    import corner as triangle

    # pull out the parameter names and flatten the thinned chains
    parnames = sample_results['theta_labels']
    
    flatchain = sample_results['chain'][:, start::thin, :]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])

    # logify mass
    if 'mass' in parnames and 'logmass' not in parnames:
        midx = np.where(np.in1d(parnames, 'mass'))[0]
        if len(midx) > 0:
            flatchain[:, midx[0]] = np.log10(flatchain[:, midx[0]])
            parnames[midx[0]] = 'logmass'
    parnames = np.array(parnames)

    # restrict to parameters you want to show
    if showpars is not None:
        ind_show = np.array([p in showpars for p in parnames], dtype=bool)
        flatchain = flatchain[:, ind_show]
        #truths = truths[ind_show]
        parnames = parnames[ind_show]

    # Make nice labels.
    niceparnames = _niceparnames(parnames)
        
    if trim_outliers:
        trim_outliers = len(parnames) * [trim_outliers]

    fig = triangle.corner(flatchain, labels=niceparnames, truths=truths,  verbose=False,
                          quantiles=[0.25, 0.5, 0.75], range=trim_outliers,
                          color='k', **kwargs)
    
    return fig
