from pylab import *
import atpy

infile='/home/share/research/luminosityfunction/testLF.fits'
lf=atpy.Table(infile)


figure()
#ax=gca()
#ax.set_xscale('log')
t=hist(lf.logL,bins=38)
clf()
bins=t[1]
yvalue=t[0]
bincenter=[]
for x in range(len(yvalue)):
    bincenter.append((bins[x]+(bins[x+1]))/2)
    #print bincenter

yerror=sqrt(t[0])    
yplot=log10(yvalue)    
plot(bincenter,yplot,'ko')
ylabel('log(N)')
xlabel('log(L)')
legend(['Sample Data Set'],loc='upper right')


yerrdown=log10(yvalue)-log10(yvalue-yerror)
yerrdown2=[]
lf2=[]
for x in range(len(yerrdown)):
    if math.isnan(yerrdown[x]):
        print 'test'
    elif yerrdown[x] > 1*10**25:
        print 'test'
    else:
        yerrdown2.append(yerrdown[x])
        lf2.append(lf[x])

yerrup=log10(yvalue+yerror)-log10(yvalue)
yerrup2=[]
lf2=[]

for x in range(len(yerrup)):
    if math.isnan(yerrup[x]):
        print 'test'
    elif yerrup[x] > 1*10**25:
        print 'test'
    elif yerrup[x]==.3010299956639812:
        print 'test'
    
    else:
        yerrup2.append(yerrup[x])
        lf2.append(lf[x])
yerrboth=zip(yerrdown,yerrup)
yerrboth=transpose(yerrboth)
errorbar(bincenter,yplot,yerrboth)

#yerror=sqrt(lf2[0])        
#plot(lf2,yplot,'ko')
        
chisqumin=10*10**5
alphabest=0.
logphistarbest=0.
loglstarbest=0.

for i in range(len(yerror)):
    if yerror[i]==0:
        yerror[i]=.1

for alpha in arange(-2,2,.1):
    for loglstar in arange(-1,1,.1):
        for logphistar in arange(-3,3,.1):
            a=log(10)*10.**logphistar
            b=(10.**array(bincenter))/(10.**loglstar)
            c=alpha+1.
            d=exp(-(10.**array(bincenter))/(10.**loglstar))
            yfit=a*(b**c)*d
            chisq=sum(log(array(yfit-yvalue)**2./array(yerror)**2.))
            if chisq<chisqumin:
                alphabest=alpha
                loglstarbest=loglstar
                logphistarbest=logphistar
                chisqumin=chisq
                yf=a*(b**c)*d
                                
plot(bincenter,log10(yf),label='Chi squared fit')
show()
print 'alpha: ' + str(alphabest)
print 'loglstar: ' + str(loglstarbest)
print 'logphistar: ' + str(logphistarbest)
print 'chisqumin: '+str(chisqumin)
