import numpy as np
import matplotlib.pyplot as plt

"""Convenience functions for reading and reconstructing results from a fitting
run, including reconstruction of the model for making posterior samples.

"""
def _niceparnames(parnames):
    """Replace parameter names with nice names."""

    old = np.array(['tau',
           'tage',
           'mass',
           'logmass',
           'logzsol',
           'dust2'])
    new = np.array([r'$\tau$ (Gyr)',
           'Age (Gyr)',
           r'$M / M_{\odot}$',
           r'$\log_{10}\,(M / M_{\odot})$',
           r'$\log_{10}\, (Z / Z_{\odot})$',
           r'$\tau_{diffuse}$'])

    niceparnames = parnames.copy().astype(new.dtype)
    for oo, nn in zip( old, new ):
        this = np.in1d(parnames, oo)
        if np.sum(this) > 0:
            niceparnames[this] = nn

    return niceparnames

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
    
    parnames = np.array(sample_results['theta_labels'], dtype='>U15')

    # logify mass
    if 'mass' in parnames:
        midx = [l == 'mass' for l in parnames]
        chain[:, :, midx] = np.log10(chain[:, :, midx])
        parnames[midx] = 'logmass'

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
    parnames = np.array(sample_results['theta_labels'], dtype='>U15')
    
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

    # Make nice labels.
    niceparnames = _niceparnames(parnames)
        
    if trim_outliers:
        trim_outliers = len(parnames) * [trim_outliers]

    fig = triangle.corner(flatchain, labels=niceparnames, truths=truths,  verbose=False,
                          quantiles=[0.25, 0.5, 0.75], range=trim_outliers, **kwargs)
    
    return fig
