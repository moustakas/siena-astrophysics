#!/usr/bin/env python
from pylab import *
import atpy
from scipy import optimize

infile = '/home/share/research/nsa/nsa_v0_1_2topcat.fits'
wise = atpy.Table(infile)
sdssflag = (wise.ISDSS >-1)
sdss2 = where(sdssflag)
sloan = wise.ISDSS[sdss2]
mag=wise.ABSMAG[sdss2]
print mag
print sloan

amag = mag[:,5]
irlum = (amag-4.83)/-2.5


t=hist(irlum, bins=50)
clf()
print t
#split t
#bins =xvalues
bins = t[1]
yvalue = t[0]

bincenters=[]
for x in range(len(yvalue)):
    bincenters.append((bins[x]+bins[x+1])/2)
    print bincenters
    print len(bincenters)
            
yerror = sqrt(t[0])
yplot = log10(yvalue)
plot(bincenters, yplot,'bo')
show()
ylabel('log(N)')
xlabel('log(L)')
yerrup=log10(yvalue+yerror)-log10(yvalue)
yerrdown=log10(yvalue)-log10(yvalue-yerror)
yerrboth=zip(yerrdown,yerrup)
yerrboth=transpose(yerrboth)
errorbar(bincenters, yplot, yerr=yerrboth)
legend(['All-Sky Survey'], loc='upper right')






    
"""alphabest=0.
logphistarbest=0.
loglstarbest=0.
chisqmin=1.*10**7

for i in range(len(yerror)):
    if yerror[i]==0:
        yerror[i]=0.1

for alpha in arange(-2,2,.1):
    for loglstar in arange(-1,1,.1):
        for logphistar in arange(-3,3,.1):
            a=log(10)*10.**logphistar
            b=(10.**array(bincenters))/(10.**loglstar)
            c=alpha+1
            d=exp(-(10**array(bincenters))/(10.**loglstar))
            yfit=a*(b**c)*d
            chisq = sum(log(array(yfit-yvalue)**2./array(yerror)**2.))
            if chisq < chisqmin:
                alphabest = alpha
                loglstarbest = loglstar
                logphistarbest = logphistar
                chisqmin=chisq
                yfit2=a*(b**c)*d

#plot(bincenters, log10(yfit2))
#show()

print 'alpha = ' + str(alphabest)
print 'loglstar = ' + str(loglstarbest)
print 'logphistar = ' + str(logphistarbest)
print 'chisqmin = ' + str(chisqmin)"""

fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
pinit = [1.0, -1.0]
out = optimize.leastsq(errfunc, pinit, args=(bincenters, yplot, yerror), full_output=1)



        


