from pylab import *
import atpy

class LF:
    def __init__(self):
        infile = '/home/share/research/luminosityfunction/testLF.fits'
        self.ldat=atpy.Table(infile)
    def plotLF(self):
        figure()
        t=hist(self.ldat.logL, bins=50)
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

LF=LF()
