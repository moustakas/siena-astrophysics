from pylab import *
import atpy

def __init__(self):
    infile='/home/share/research/luminosityfunction/testLF.fits'
    self.ldat=atpy.Table(infile)
    
def plotlum(self):
     figure()
     t=hist(self.ldat.logL,bins=100)
     clf()
     bins=t[1]
     yvalue=t[0]
     yerror=sqrt(t[0])
     bincenter=[]
     for x in range(len(yvalue)):
         bincenter.append((bins[x]+(bins[x+1]))/2)
         print bincenter
     yplot=log10(yvalue)    
     plot(bincenter,yplot,'ko')
     yerrup=log10(yvalue+yerror)-log10(yvalue)
     yerrdown=log10(yvalue)-log10(yvalue-yerror)
     yerrboth=zip(yerrdown,yerrup)
     yerrboth=transpose(yerrboth)
     errorbar(bincenter,yplot,yerrboth)
     #lummy=arange(-.01,5.,.1)
     Lvalues=linspace(-2.,1.,100)
     Llin=10**(Lvalues)
     schfun=(Llin**(-.25))*exp(-Llin)*55
     plot(Lvalues,log10(schfun))
test=test()
