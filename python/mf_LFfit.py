#!/usr/bin/env python
from pylab import *
import atpy
import numpy
 


infile = '/home/share/research/luminosityfunction/testLF.fits'
lf=atpy.Table(infile)

figure()
t=hist(lf.logL, bins=50)
clf()
print t
#split t
#bins =xvalues
bins = t[1]
yvalue = t[0]

bincenters=[]
for x in range(len(yvalue)):
    bincenters.append((bins[x]+bins[x+1])/2)

#print bincenters
#print len(bincenters)            
yerror =(sqrt(t[0]))
yplot = (log10(yvalue))
plot(bincenters, yplot,'bo')
ylabel('log(N)')
xlabel('log(L)')


"""yerrdown2=[]
lf2=[]
bincenters2=[]
yerrup2=[]
for x in range(len(yerrdown)):
    if math.isnan(yerrdown[x]):
        print 'hey'
    elif yerrdown[x] > 1*10**25:
        print 'hey'
    else:
        yerrdown2.append(yerrdown[x])
        lf2.append(lf[x])
        bincenters2.append(lf[x])


for x in range(len(yerrup)):
    if math.isnan(yerrup[x]):
        print 'text1'
    elif yerrup[x] > 1*10**25:
        print 'text2'
    elif yerrup[x]== .3010299956639812:
        print 'text3'
    else:
        yerrup2.append(yerrup[x])
        lf2.append(lf[x])
        bincenters2.append(lf[x])

lf2=hist(lf.logL, bins=46)
clf()
bins2=lf2[1]
yvalue2=lf2[0]
bincenters2=[]
for x in range(len(yvalue2)):
    bincenters2.append((bins2[x]+bins2[x+1])/2)
yerror2=(sqrt(lf2[0]))
yplot2=(log10(yvalue2))
plot(bincenters2, yplot2, 'r*')"""


yerrdown=log10(yvalue)-log10(yvalue-yerror)
yerrup=log10(yvalue+yerror)-log10(yvalue)
yerrboth2=zip(yerrup2,yerrdown2)
yerrboth2=transpose(yerrboth2)
errorbar(bincenters2, yplot2, yerr=yerrboth2)



#ax=gca()
#ax.set_yscale('log')
#ax=gca()
#ax.set_xscale('log')
#Lvalues = linspace(-2, 1, 100)
#Llinear = 10**(Lvalues)
#schect = 55*(Llinear**(-0.25))*exp(-Llinear)
#plot(Lvalues, log10(schect))






alphabest=0.
logphistarbest=0.
loglstarbest=0.
chisqmin=1.*10**7

for i in range(len(yerror2)):
    if yerror2[i]==0:
        yerror2[i]=0.1

for alpha in arange(-2,2,.1):
    for loglstar in arange(-1,1,.1):
        for logphistar in arange(-3,3,.1):
            yfit = (log(10)*10**(logphistar))*((yvalue/10**loglstar)**alpha+1)*exp(yvalue/10**loglstar)
            chisq = sum((yfit-yvalue)**(2))/(yerror**2)
            if chisq < chisqmin:
                alphabest = alpha
                loglstarbest = loglstar
                logphistarbest = logphistar
                print alpha
                yfit2 = (log(10)*10**(logphistar))*((yvalue/10**loglstar)**alpha+1)*exp(yvalue/10**loglstar)
plot(bincenters, yfit2)
show()































































































			
        
