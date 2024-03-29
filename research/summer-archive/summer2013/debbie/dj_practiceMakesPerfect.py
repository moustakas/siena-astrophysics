

#! /user/bin/env python
from pylab import *
import atpy
from LCScommon import *



class cluster:
    def __init__(self,clustername):
        infile='/home/share/research/LocalClusters/NSAmastertables/'+clustername+'_NSAmastertable.fits'
        self.mdat=atpy.Table(infile)
        self.clusterRA=clusterRA[clustername]

    def plotradec(self):
        figure()
        plot(self.mdat.RA,self.mdat.DEC,'b.')
    def plotrs(self):
        figure()
        hist(self.mdat.ZDIST,bins=20,histtype='step')
        
    def plotlf(self):
        #figure()
        hist(self.mdat.ABSMAG[:,3],bins=100)

    def plotvel(self):
        dr =sqrt((self.mdat.RA-self.clusterRA)**2+(self.mdat.DEC)**2)
        ax=gca()
        ax.set_yscale('log')
        plot(dr,self.mdat.VDISP)

    def lum(self,myaxes,myaxes_dummy):
        #figure()
        #self.absmag=self.mdat.ABSMAG[:,3]
        #self.lum=(self.absmag-4.83)/-2.5
        #ylabel=('Number of Galaxies')
        #xlabel=('log(L/Lsun)')
        #hist(self.lum,bins=100)       
        
        
        t=myaxes_dummy.hist(self.mdat.ABSMAG[:,3],bins=100)
        #clf()
        bins=t[1]
        yvalue=t[0]
        yerror=sqrt(t[0])
        bincenter=[]
        for x in range(len(yvalue)):
            bincenter.append((bins[x]+(bins[x+1]))/2)
            #print bincenter
        #apple=gca()    
        #apple.set_xvalues('log')
        #myaxes.set_xscale('log')
        yplot=log10(yvalue)    
        myaxes.plot(bincenter,yplot,'ko')
        yerrup=log10(yvalue+yerror)-log10(yvalue)
        yerrdown=log10(yvalue)-log10(yvalue-yerror)
        yerrboth=zip(yerrdown,yerrup)
        yerrboth=transpose(yerrboth)
        myaxes.errorbar(bincenter,yplot,yerrboth,axes=myaxes)
    

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
    fig = figure(1)
    clf()
    fig_dummy = figure(2)
    clf()
    ax_dummy = fig_dummy.add_subplot(1,1,1)
    while f<10:
        ax = fig.add_subplot(3,3,f)
        xlabel('g band log(L/L$_\odot$)')
        ylabel('Number of Galaxies')
        (clusters[f-1]).lum(ax,ax_dummy)
        ax.legend([clusternames[f-1]])
        #lummy=arange(-.01,5.,.1)
        #schfun=6.*(((lummy)**(-.25+1))*(exp(-lummy)))
        #ax.plot(lummy,schfun)
        #lax=gca()
        #lax.set_xscale('log')
        Lvalues=linspace(-2.,1.,100)
        Llin=10**(Lvalues)
        #schfun=(Llin**(-.25))*exp(-Llin)*55
        #ax.plot(Lvalues,log10(schfun))
        f=f+1
    #figure(1)
    
    ax=fig.gca()
    figure(1)
    ax2=gca()
    print ax,ax2
    text(-.5,-.4,'g band log(L/L$_\odot$)',transform=ax.transAxes,horizontalalignment='center',fontsize=20)
    text(-2.8,1.5,'log Number of Galaxies',transform=ax.transAxes,verticalalignment='center',rotation=90,fontsize=20)
def plotboth():
    f=1
    while f<10:
        subplot(3,3,f)
        xlabel('Mg')
        ylabel('Number of Galaxies')
        hist(clusters[f-1].mdat.ABSMAG[:,3],bins=20)
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
        xlabel('u"\u03b8"')
        plot(lummy,schfun)
        ax=gca()
        ax.set_yscale('log')
        ax=gca()
        ax.set_xscale('log')
       # schfun=(8.*10**-3)((1/((3.9*10**33)(1.4*10**10))**alpha)*exp(-1/((3.9*10**33)(1.4*10**10)))
        #plot(,schfun,'r')
        f=f+1

