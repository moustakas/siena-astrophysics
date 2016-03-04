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
def lnlike(theta,x,y,yerr):
	 m, b, lnf = theta
	 model = m * x + b
	 inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
	 return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def lnprior(theta):
	m, b, lnf = theta
	print m	
	if type(m) == np.ndarray:
		print np.shape(m)
		import pdb ; pdb.set_trace()
	if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
		return 0.0
	return -np.inf

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)


#The "true" parameters.
m_true = -0.9594
b_true = 4.294
f_true = 0.534

# Generate some synthetic data from the model.
N = 50
x = np.sort(10*np.random.rand(N))
yerr = 0.1+0.5*np.random.rand(N)
y = m_true*x+b_true
y += np.abs(f_true*y) * np.random.randn(N)
y += yerr * np.random.randn(N)

nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))
m_ml, b_ml, lnf_ml = result["x"]


nwalkers, ndim = 50, 3
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob, args = (x,y,yerr))
sampler.run_mcmc(pos,500)
#import pdb
#pdb.set_trace()

samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
fig = corner.corner(samples, labels=["$m$", "$b$", "$\ln\,f$"],
                      truths=[m_true, b_true, np.log(f_true)])
fig.savefig("triangle.png")

plt.show()

