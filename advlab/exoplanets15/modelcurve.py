import matplotlib.pyplot as plt
import numpy as np
import batman
import matplotlib
import kplr
################Model Light Curve##############################

def lnlike(mflux,rflux,rferr):
	 return -.5*np.sum((rflux-mflux)**2/rferr**2)/2
	
time, flux, ferr = np.loadtxt("k30-c.txt",unpack = True)

client = kplr.API()
koi = client.koi("806.02")

Rearth = 6.371e6
Rsun = 6.957e8
autosun = Rsun / 1.5e11

# Set the parameters for the model curve
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
m = batman.TransitModel(params, t/24.0)

modelcurve = m.light_curve(params)
print lnlike(modelcurve,flux,ferr)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(time,flux,'gs')
y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.yaxis.set_major_formatter(y_formatter)
ax.plot(t,modelcurve,'-r')
plt.ylim(np.min(flux) - .01 ,np.max(flux))
#plt.xlim(-0.1,0.1)
plt.show()
