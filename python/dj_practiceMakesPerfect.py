

#! /user/bin/env python
from pylab import *
import atpy



class cluster:
    def __init__(self,clustername):
        infile='/home/share/research/LocalClusters/NSAmastertables/'+clustername+'_NSAmastertable.fits'
        self.mdat=atpy.Table(infile)
   

    def plotradec(self):
        figure()
        plot(self.mdat.RA,self.mdat.DEC,'b.')
    def plotrs(self):
        figure()
        hist(self.mdat.ZDIST,bins=20,histtype='step')
        
    def plotlf(self):
        #figure()
        hist(self.mdat.ABSMAG[:,3],bins=100)

    def lum(self):
        #figure()
        self.absmag=self.mdat.ABSMAG[:,3]
        self.lum=(self.absmag-4.83)/-2.5
        ylabel=('Number of Galaxies')
        xlabel=('log(L/Lsun)')
        hist(self.lum,bins=100)       
    
    

mkw11=cluster('MKW11')
mkw8=cluster('MKW8')
coma=cluster('Coma')
herc=cluster('Hercules')
ngc=cluster('NGC6107')
awm4=cluster('AWM4')
a1367=cluster('A1367')
a2052=cluster('A2052')
a2063=cluster('A2063')
clusters=[mkw11,mkw8,coma,herc,ngc,awm4,a1367,a2052,a2063]
clusternames=['MKW11','MKW8','Coma','Hercules','NGC6107','AWM4','A1367','A2052','A2063']

def lumall():
    f=1
    while f<10:
        subplot(3,3,f)
        xlabel('g band log(L/L$_\odot$')
        ylabel('Number of Galaxies')
        (clusters[f-1]).lum()
        legend([clusternames[f-1]])
        f=f+1
    
def plotboth():
    f=1
    while f<10:
       # subplot(3,3,f)
        xlabel('Mg')
        ylabel('Number of Galaxies')
        hist(clusters[f-1].mdat.ABSMAG[:,3],bins=20,histtype='step')
        #ax=gca()
        #ax.set_yscale('log')
        legend([clusternames[f-1]])
        f=f+1
    
def fun():
    f=0
    alpha=[-.5,-.75,-1.]
    lummy=arange(.01,5.,.1)
    figure()
    while f<3:
        schfun=(8.0*10.**-3.)*((lummy)**(alpha[f]))*(exp(-lummy))
        
        plot(lummy,schfun)
        ax=gca()
        ax.set_yscale('log')
        ax=gca()
        ax.set_xscale('log')
       # schfun=(8.*10**-3)((1/((3.9*10**33)(1.4*10**10))**alpha)*exp(-1/((3.9*10**33)(1.4*10**10)))
        #plot(,schfun,'r')
        f=f+1
