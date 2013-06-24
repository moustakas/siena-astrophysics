from pylab import *
import atpy

class LF:
    def __init__(self):
        infile = '/home/share/research/luminosityfunction/testLF.fits'
        self.ldat=atpy.Table(infile)
    def plotLF(self):
        figure()
        t=hist(self.ldat.logL, bins=100)
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
        ylabel('log(N)')
        xlabel('log(L)')
        yerrup=log10(yvalue+yerror)-log10(yvalue)
        yerrdown=log10(yvalue)-log10(yvalue-yerror)
        yerrboth=zip(yerrdown,yerrup)
        yerrboth=transpose(yerrboth)
        errorbar(bincenters, yplot, yerr=yerrboth)
        #ax=gca()
        #ax.set_yscale('log')
        #ax=gca()
        #ax.set_xscale('log')
        Lvalues = linspace(-2, 1, 100)
        Llinear = 10**(Lvalues)
        schect = 55*(Llinear**(-0.25))*exp(-Llinear)
        plot(Lvalues, log10(schect))
       
        
        
        

LF=LF()
