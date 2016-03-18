#!/usr/bin/python

import sys
import kplr
import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve, Box1DKernel
import batman
import copy
import corner
import emcee

time = []
flux = []
ferr = []
quality = []

client = kplr.API()
koi = client.koi(806.02)

######Functions for Optimization##########################

def lnlike(theta,params,time,flux,ferr):
    myparams = copy.copy(params)
    #print theta
    myparams.per = theta[0]
    myparams.rp = theta[1]
    m = batman.TransitModel(myparams, t/24.0)
    mflux = m.light_curve(myparams)
    return -0.5*np.sum((flux-mflux)**2/ferr**2)

def lnprior(theta):

    per = theta[0]
    rp = theta[1]

    if ((per>50.00)*(per<70.00)) and ((rp>0.1)*(rp<0.5)):
        return 0.0
   
    return -np.inf

def lnprob(theta,params,time,flux ,ferr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta,params,time,flux,ferr)

####################Get Light Curves####################

lcs = koi.get_light_curves(short_cadence = False)

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
print flux
plt.figure()
plt.plot(time,flux,'go')

###############Normalized Stacked Light Curve########################

midcurve = koi.koi_time0bk-67 
per = koi.koi_period 
dur = koi.koi_duration/24
neclipse = 7 
ndeg = 2

stacktime = []
normalflux = []
stacknorm = []
stackferr = []

for ii in range(neclipse):
	wlong = np.where( (time > (midcurve + ii*per - dur - 0.75)) * 
		(time < (midcurve + ii *per + dur + 0.75)) * 1)

	thistime = time[wlong]
	thisferr = ferr[wlong]
	thisflux = flux[wlong]

	thisflux[np.where(np.isnan(thisflux))] = np.median(thisflux[~np.isnan(thisflux)]) 
	thisferr[np.where(np.isnan(thisferr))] = np.median(thisferr[~np.isnan(thisferr)]) 

	stacktime.append((thistime - ii * per - midcurve)*24.0)
	
	wshort = np.where( (thistime > (midcurve + ii*per - 1.3* dur)) *
                (thistime < (midcurve + ii * per + 1.3*dur)) * 1)	
	ivar = 1 / thisferr**2
	ivar[wshort] = 0

	coeff = np.polyfit(thistime, thisflux, ndeg,  w = ivar)

#	print("thistime", thistime)
#	print("thisflux", thisflux)
#	print("ivar", ivar)
#	print("coeff", coeff)
	
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
#############################################################

stacknorm = np.concatenate(stacknorm)
stacktime = np.concatenate(stacktime)
stackferr = np.concatenate(stackferr)
np.savetxt("k30-c.txt",zip(stacktime,stacknorm, stackferr))

################Model Light Curve##############################

Rearth = 6.371e6
Rsun = 6.955e8

true_per =  koi.koi_period
true_rp = koi.koi_prad * (Rearth / Rsun) / koi.koi_srad
params = batman.TransitParams()
params.t0 = 0.
params.per = true_per
params.rp = true_rp
params.a = koi.koi_sma * 1.496E11 / Rsun
params.inc = koi.koi_incl
params.ecc = koi.koi_eccen
params.w = 90 #koi.koi_longp
params.limb_dark = 'nonlinear'
params.u = [koi.koi_ldm_coeff1,koi.koi_ldm_coeff2,koi.koi_ldm_coeff3,koi.koi_ldm_coeff4]

t = np.linspace(-20,20,len(stacknorm))

###################Optimization######################

nwalkers, ndim = 100, 2
pos = np.zeros((nwalkers,ndim))

pos[:,0] = np.random.uniform(3.0, 6.0, nwalkers)
pos[:,1] = np.random.uniform(0, 0.5, nwalkers)

#print(pos)

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(params,stacktime,stacknorm,stackferr), threads = 3)
sampler.run_mcmc(pos,2000)

samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
np.savetxt('k30c-samples.txt',zip(samples[0],samples[1]))

fig = corner.corner(samples,labels = ["per","a"], truths = [true_per,true_rp])
plt.show()

####################total plot########################
plt.figure()
#plt.plot(t,modelcurve)
plt.plot(stacktime, stacknorm,'go')
plt.ylabel('Normalized Flux')
plt.xlabel('Time (Hours)')
plt.title('Stacked and Normalized Curve')
plt.xlim((-9.0, 9.0))

plt.show()

