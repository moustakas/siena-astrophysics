import matplotlib.pyplot as plt
import numpy as np
import batman
import matplotlib
import kplr
import emcee
import scipy.optimize as op
import sys
import corner
import copy


###################Set up functions to calculate probablilities##################
def lnlike(theta,params,time,flux,ferr):
    params.per = theta
    m = batman.TransitModel(params, t/24.0)
    mflux = m.light_curve(params)
    plt.plot(time,mflux, c = 'b')
    return -0.5*np.sum((flux-mflux)**2/ferr**2)

def lnprior(theta):
    per = theta
    print per
	
    if type(per) == np.ndarray:
        import pdb ; pdb.set_trace()
	
    if (per>1.00)*(per<10.00):
        return 0.0
    
    return -np.inf

def lnprob(theta,params,time,flux ,ferr):
    #print 'theta', theta
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta,params,time,flux,ferr)

time, flux, ferr = np.loadtxt("k30-c.txt",unpack = True)

client = kplr.API()
koi = client.koi("806.02")

Rearth = 6.371e6
Rsun = 6.957e8
autosun = Rsun / 1.5e11

################ Set parameters for the model light curve####################
params = batman.TransitParams()
params.t0 = 0.
params.per = koi.koi_period
params.rp = koi.koi_prad * (Rearth / Rsun) / koi.koi_srad
params.a = koi.koi_sma * 1.496E11 / Rsun 
params.inc = koi.koi_incl
params.ecc = koi.koi_eccen
params.w = 80#koi.koi_longp
params.limb_dark = "nonlinear"
params.u = [koi.koi_ldm_coeff1,koi.koi_ldm_coeff2,koi.koi_ldm_coeff3,koi.koi_ldm_coeff4]
#print(params.a)

t = np.linspace(-20,20,len(flux))

nwalkers, ndim = 30, 1
pos = np.random.uniform(1.0,9.0,nwalkers) 
sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob, args = (params,time,flux,ferr))
sampler.run_mcmc(pos,500)
#import pdb
#pdb.set_trace()

#samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
plt.plot(samples)

plt.show()


'''
period = np.linspace(1.5,9.5,20)
prob = np.zeros_like(period)
myparams = params

plt.plot(time,flux,'gs')
for ii,per in enumerate(period):
    myparams.per = per
    prob[ii] = lnprob(myparams,time,flux,ferr)
    print(per,prob[ii])
plt.show()

plt.plot(period,prob,marker = 's',ls = '-')
plt.show()

sys.exit(1)

m = batman.TransitModel(params, t/24.0)
modelcurve = m.light_curve(params)

###Attempt to use emcee################
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, ["prad"], args=(modelcurve, flux, ferr))
m_ml, b_ml, lnf_ml = result["prad"]

ndim, nwalkers = 3, 100
pos = [result["prad"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob, args = (modelcurve,flux,ferr))
sampler.run_mcmc(pos,500)
########################################

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(time,flux,'gs')
y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.yaxis.set_major_formatter(y_formatter)
ax.plot(t,modelcurve,'-r')
plt.ylim(np.min(flux) - .01 ,np.max(flux))
#plt.xlim(-0.1,0.1)
plt.show()
'''
