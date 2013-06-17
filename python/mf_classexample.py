#!/usr/bin/env python
from pylab import *
import atpy


class cluster:
    def __init__(self,clustername):
        infile='/home/rfinn/research/LocalClusters/NSAmastertables/'+clustername+'_NSAmastertable.fits'
        self.mdat=atpy.Table(infile)
    def plotradec(self):
        figure()
        plot(self.mdat.RA,self.mdat.DEC,'b.')
    def plotredshift(self):
        figure()
        hist(self.mdat.ZDIST)
        title('RedShift MKW11')
        

class nasasloan:
    def __init__(self):
        infile='/home/share/research/nsa/nsa_v0_1_2topcat.fits'
        #infile='/home/share/research/nsa/NSA_LCSregion.fits'
        self.ndat=atpy.Table(infile)
    def gbandabsmag(self):
        figure()
        hist(self.ndat.ABSMAG[:,3],bins=100)
        #title('Absolute Magnitude( g-band)from WISE ') 
        ax=gca()
        ax.set_yscale('log')
        xlabel('M')
        ylabel('N')
    def gbandlum(self): 
        figure()
        self.amag=self.ndat.ABSMAG[:,3]
        print self.amag
        self.luminosity=(self.amag-4.83)/-2.5
        print 'this is a test '
        print self.luminosity
        hist(self.luminosity,bins=500,range=[5,12])
        #hist(self.luminosity,bins=100,)
        # ax=gca()
        # ax.set_yscale('log')
        xlim(5,12)
        legend(['WISE Data'], loc='upper left')
        xlabel(' log(g-band L/L$_\odot$)') 
        ylabel('N')
    def plotLF(self):
        self.amag=self.ndat.ABSMAG[:,3]
        self.luminosity=(self.amag-4.83)/-2.5
        figure()
        t=hist(self.luminosity, bins=50)
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
        legend(['All-Sky Survey'], loc='upper right')



nsa=nasasloan()
    
