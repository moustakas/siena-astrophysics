#!/usr/bin/python

import sys
import kplr
import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve, Box1DKernel
import batman

time = []
flux = []
ferr = []
quality = []

client = kplr.API()
koi = client.koi("806.02")

################Model Light Curve##############################

Rearth = 6.371e6
Rsun = 6.955e8
autosun = Rsun / 1.5e11

# Set the parameters for the model curve
params = batman.TransitParams()
params.t0 = 0. 
params.per = koi.koi_period
params.rp = koi.koi_prad * (Rearth / Rsun) / koi.koi_srad
params.a = koi.koi_sma * (autosun / koi.koi_srad)
params.inc = koi.koi_incl
params.ecc = koi.koi_eccen
params.w = 90#koi.koi_longp
params.limb_dark = "nonlinear"
params.u = [koi.koi_ldm_coeff1,koi.koi_ldm_coeff2,koi.koi_ldm_coeff3,koi.koi_ldm_coeff4]

t = np.linspace(-0.2392/2,0.2392/2,678)
m = batman.TransitModel(params, t)
modelcurve = m.light_curve(params) 

#plt.figure()
#plt.plot(t,modelcurve)
#plt.ylim(.890,1)

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

midcurve = 176.904 
per = 60.3251 
dur = 0.2392
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

	print("thistime", thistime)
	print("thisflux", thisflux)
	print("ivar", ivar)
	print("coeff", coeff)
	
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

plt.figure()
stacknorm = np.concatenate(stacknorm)
stacktime = np.concatenate(stacktime)
stackferr = np.concatenate(stackferr)
np.savetxt("k30-c.txt",zip(stacktime,stacknorm, stackferr))
plt.plot(t,modelcurve)
plt.plot(stacktime, stacknorm,'go')
plt.ylabel('Normalized Flux')
plt.xlabel('Time (Hours)')
plt.title('Stacked and Normalized Curve')
plt.xlim((-9.0, 9.0))

plt.show()

