#! /usr/bin/env python

import os
import numpy as np
from astropy.io import fits
import seaborn as sns
import matplotlib.pyplot as plt

def medxbin(x,y,binsize,minpts=20,xmin=None,xmax=None):
    """
    Compute the median (and other statistics) in fixed bins along the x-axis. 
    """
    import numpy as np
    from scipy import ptp

    # Need an exception if there are fewer than three arguments.

    if xmin==None:
        xmin = x.min()
    if xmax==None:
        xmax = x.max()
    #print(xmin,xmax)

    nbin = long(ptp(x)/binsize)
    bins = np.linspace(xmin,xmax,nbin)
    idx  = np.digitize(x,bins)
    #print(nbin, bins, xmin, xmax)

    stats = np.zeros(nbin,[('npts','f4'),('median','f8'),('sigma','f8'),('iqr','f8')])
    for kk in np.arange(nbin):
        npts = len(y[idx==kk])
        if npts>minpts:
            stats['npts'][kk] = npts
            stats['median'][kk] = np.median(y[idx==kk])
            stats['sigma'][kk] = np.std(y[idx==kk])
            stats['iqr'][kk] = np.subtract(*np.percentile(y[idx==kk],[75, 25]))

    # Remove bins with too few points.
    good = np.nonzero(stats['median'])
    stats = stats[good]

    return bins[good], stats

def main():
    """This script simply takes input cluster richness values and outputs halo mass.

    Also outputs a figure of M(BCG) vs M(Halo).

"""
    in_dir = os.getenv('HOME')+'/redmapper/'
    
    isedfit = fits.getdata(in_dir+'decals-isedfit.fits')
    redmapper = fits.getdata(in_dir+'decals-redmapper.fits')
    
    redshift = isedfit['z']
    mbcg = isedfit['mstar_avg']
    lam = redmapper['lambda_chisq']
    mhalo = np.log10((10.0**14.0)*np.exp(1.72+1.08*np.log(lam/60.0)))
    midz = np.where(((redshift>0.2)*1) & ((redshift<0.4)*1))



    
    #M(halo) vs M(BCG)
    sns.set(style='white',font_scale=1.5)
    fig = plt.figure(figsize=(8,6))
    ax = fig.gca()
    sns.kdeplot(mhalo[midz], mbcg[midz], shade=True)
    sns.despine(fig=fig)
    #plt.scatter(mhalo, mbcg, c='b')
    plt.xlabel('$log_{10}\ (M_{halo}/M_{\odot})$')
    plt.ylim(10.5,12.5)
    plt.xlim(13.5,14.8)
    plt.ylabel('$log_{10}\ (M_{BCG}/M_{\odot})$')

    m, b = np.polyfit(mhalo[midz], mbcg[midz], 1)
    alpha = str(b)[0:4]
    beta = str(m)[0:5]
    xmodel = np.linspace(13.5, 14.8, 50)
    ymodel = np.polyval(np.polyfit(mhalo[midz], mbcg[midz], 1), xmodel)
    plt.plot(xmodel, ymodel, 'k-')
    fig.subplots_adjust(bottom=0.2, left=0.15)
    ax.text(0.95, 0.90, '$M_{BCG}\ \propto \, (M_{halo})^{'+beta+'}$',
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=14)
    ax.text(0.95, 0.10, '0.2<Redshift<0.4', horizontalalignment='right',
            verticalalignment='bottom', transform=ax.transAxes,
            fontsize=14)
    plt.savefig(in_dir+'mbcg_mhalo.jpg',clobber=True)
    plt.close()




    #Galaxy Formation Efficiency
    sns.set(style='white',font_scale=1.5)
    fig2 = plt.figure(figsize=(8,6))
    
    plt.scatter(mhalo, mbcg-mhalo, c='b')
    sns.despine(right=True, top=True)
    plt.xlabel('$log_{10}\ (M_{halo}/M_{\odot})$')
    plt.ylabel('$log_{10}\ (M_{BCG}/M_{halo})$')
    fig2.subplots_adjust(bottom=0.2, left=0.15)
    plt.savefig(in_dir+'galaxy_efficiency.jpg',clobber=True)
    plt.close()



    
    #M(BCG) vs redshift
    sns.set_style('white')
    lo = np.where((mhalo<14.0)*1)
    mid = np.where(((mhalo>14.0)*1) & ((mhalo<14.5)*1))
    hi = np.where(mhalo>14.5)
    binsz = 0.04
    lobins, lostats = medxbin(redshift[lo], mbcg[lo], binsz)
    midbins, midstats = medxbin(redshift[mid], mbcg[mid], binsz)
    hibins, histats = medxbin(redshift[hi], mbcg[hi], binsz)

    fig3 = plt.figure(figsize=(8,6))
    #sns.kdeplot(redshift[lo], mbcg[lo], cmap="Blues", legend=True)
    #sns.kdeplot(redshift[mid],mbcg[mid], cmap="Reds", legend=True)
    #sns.kdeplot(redshift[hi],mbcg[hi], cmap="Greens_r", legend=True)
    
    #plt.plot(lobins, lostats['median'])
    #plt.plot(midbins, midstats['median'])
    #plt.plot(hibins, histats['median'])
    plt.errorbar(lobins, lostats['median'], yerr=lostats['sigma']/np.sqrt(lostats['npts']), c='b')
    plt.errorbar(midbins, midstats['median'], yerr=midstats['sigma']/np.sqrt(midstats['npts']), c='r')
    plt.errorbar(hibins, histats['median'], yerr=histats['sigma']/np.sqrt(histats['npts']), c='g')
    plt.xlabel('Redshift')
    plt.ylabel('$log_{10}\ (M_{BCG}/M_{\odot})$')
    plt.ylim(10.5,12.5)
    plt.xlim(0.05,0.73)
    fig3.subplots_adjust(bottom=0.2, left=0.15)
    plt.legend(['low halo', 'mid halo', 'high halo'])
    plt.savefig(in_dir+'3_halo.jpg',clobber=True)

    
if __name__ == '__main__':
    main()
