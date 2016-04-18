#!/usr/bin/python

import sys
import kplr
import matplotlib.pyplot as plt
import numpy as np
import batman
import copy
import corner
import emcee
import argparse

Rearth = 6.371e6
Rsun = 6.957e8

def get_lightcurves(koi='806.02'):
    time = []
    flux = []
    ferr = []
    quality = []
    
    client = kplr.API()
    koidat = client.koi(koi)

    lcs = koidat.get_light_curves(short_cadence = False)

    # download all of the different light curves
    for lc in lcs:
    	with lc.open() as f:
		hdu_data = f[1].data
		time.append(hdu_data["time"])
		flux.append(hdu_data["sap_flux"])
		ferr.append(hdu_data["sap_flux_err"])
		quality.append(hdu_data["sap_quality"])

def normalize(koidat):
    time = []
    flux = []
    ferr = []
    quality = []
    
    lcs = koidat.get_light_curves(short_cadence = False)
    name = koidat.kepler_name

    # download all of the different light curves
    for lc in lcs:
        with lc.open() as f:
                hdu_data = f[1].data
                time.append(hdu_data["time"])
                flux.append(hdu_data["sap_flux"])
                ferr.append(hdu_data["sap_flux_err"])
                quality.append(hdu_data["sap_quality"])

    time = np.concatenate(time)
    time = time - 67

    #put all of the different eclipses together
    flux = np.concatenate(flux)
    ferr = np.concatenate(ferr)
    quality = np.concatenate(quality)
   
    plt.figure()
    plt.plot(time,flux,'go')

    midcurve = koidat.koi_time0bk-67 #subtract 67 to match the time scale
    per = koidat.koi_period 
    dur = koidat.koi_duration/24   #converts the duration into days
    neclipse = 7 
    ndeg = 2

    stacktime=[]
    normalflux=[]
    stacknorm=[]
    stackferr=[]

    for ii in range(neclipse):
        wlong = np.where( (time > (midcurve + ii*per - dur - 0.75)) *
                (time < (midcurve + ii *per + dur + 0.75)) * 1)

        thistime=time[wlong]
        thisferr=ferr[wlong]
        thisflux=flux[wlong]

        thisflux[np.where(np.isnan(thisflux))] = np.median(thisflux[~np.isnan(thisflux)])
        thisferr[np.where(np.isnan(thisferr))] = np.median(thisferr[~np.isnan(thisferr)])

        temptime = (thistime - ii * per - midcurve)*24.0

        stacktime.append((thistime - ii * per - midcurve)*24.0)

        wshort = np.where( (thistime > (midcurve + ii*per - 1.3* dur)) *
                (thistime < (midcurve + ii * per + 1.3*dur)) * 1)
        
        ivar = 1 / thisferr**2
        ivar[wshort] = 0

        coeff = np.polyfit(thistime, thisflux, ndeg,  w = ivar)

        fit = np.polyval(coeff,thistime)

        normalflux = thisflux / fit
        normalferr = thisferr / fit

        stacknorm.append(normalflux)
        stackferr.append(normalferr)

        plt.figure()
        plt.errorbar(thistime,thisflux,thisferr)
        plt.plot(thistime, fit, "g-")
        plt.plot(thistime[wshort], thisflux[wshort], "ro")
        plt.ylabel('Flux')
        plt.xlabel('Time(Days)')
        plt.title('Flux with Line of Best Fit (Light Curve Masked)')

    stacknorm = np.concatenate(stacknorm)
    stacktime = np.concatenate(stacktime)
    stackferr = np.concatenate(stackferr)
    np.savetxt(name + 'norm',zip(stacktime,stacknorm, stackferr))

def lnlike(theta,params,time,flux,ferr):
    myparams = copy.copy(params)
    myparams.per = theta[0]
    myparams.rp = theta[1]
    t = np.linspace(-20,20,len(flux))
    m = batman.TransitModel(myparams, t/24.0)
    mflux = m.light_curve(myparams)
    return -0.5*np.sum((flux-mflux)**2/ferr**2)

def lnprior(theta):
    per = theta[0]
    rp = theta[1]

    if ((per>50.00)*(per<70.00)) and ((rp>0.05)*(rp<0.25)):
        return 0.0
   
    return -np.inf

def lnprob(theta,params,time,flux ,ferr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta,params,time,flux,ferr)

def optimize(koidat):
    name = koidat.kepler_name
    stacktime, stacknorm, stackferr = np.loadtxt(name + 'norm',unpack = True)

    true_per = koidat.koi_period
    true_rp = koidat.koi_prad * (Rearth / Rsun) 

    params = batman.TransitParams()
    params.t0 = 0.
    params.per = true_per
    params.rp = true_rp
    params.a = koidat.koi_sma * 1.496E11 / Rsun
    params.inc = koidat.koi_incl
    params.ecc = koidat.koi_eccen
    params.w = 90 #koidat.koi_longp
    params.limb_dark = 'nonlinear'
    params.u = [koidat.koi_ldm_coeff1,koidat.koi_ldm_coeff2,koidat.koi_ldm_coeff3,koidat.koi_ldm_coeff4]

    t = np.linspace(-20,20,len(stacknorm))
    
    nwalkers, ndim = 50, 2

    pos = np.zeros((nwalkers,ndim))
    pos[:,0] = np.random.uniform(40.0, 80.0, nwalkers)
    pos[:,1] = np.random.uniform(0.0, 0.3, nwalkers)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(params,stacktime,stacknorm,stackferr), threads = 3)
    sampler.run_mcmc(pos,250)

    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    np.savetxt('samples.txt',samples)

##########Creation of Final Plots######################
def finalplots(samples,koidat):
    sample = np.loadtxt(samples,unpack = False)  
    true_per = koidat.koi_period
    true_rp = koidat.koi_prad * (Rearth / Rsun) 
    fig = corner.corner(sample,labels = ["per","rp"], truths = [true_per,true_rp])

    params = batman.TransitParams()
    params.t0 = 0.
    params.per = true_per
    params.rp = true_rp
    params.a = koidat.koi_sma * 1.496E11 / Rsun
    params.inc = koidat.koi_incl
    params.ecc = koidat.koi_eccen
    params.w = 90 #koidat.koi_longp
    params.limb_dark = 'nonlinear'
    params.u = [koidat.koi_ldm_coeff1,koidat.koi_ldm_coeff2,koidat.koi_ldm_coeff3,koidat.koi_ldm_coeff4]

    plt.show()

####################################################

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-k','--koi', type=str, default='806.02', help='koi')
    parser.add_argument('--get-lightcurves',action='store_true', help='get light curves')
    parser.add_argument('--normalize', action = 'store_true', help='Normalize the light curves')
    parser.add_argument('--optimize', action='store_true',help='Fit with MCMC')
    parser.add_argument('--final_plot', action='store_true', help='Produce triangle/final plots')

    args = parser.parse_args()
    if args.koi is None:
        parser.print_help()
        sys.exit(1)
    
    koidat = kplr.API().koi(args.koi)
    
    if args.normalize:
        normalize(koidat)
    if args.optimize:
        optimize(koidat)
    if args.final_plot:
        finalplots('samples.txt',koidat) 

if __name__ == '__main__':
    main()

