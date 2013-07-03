#! /user/bin/env python
from pylab import *
import atpy
from LCScommon import *

class cluster:
    def __init__(self,clustername):
        infile='/home/share/research/LocalClusters/NSAmastertables/'+clustername+'_NSAmastertable.fits'
        self.mdat=atpy.Table(infile)
        self.clusterRA=clusterRA[clustername]
        self.clusterDec=clusterDec[clustername]
        self.clustervel=clustervel[clustername]

    def plotvel(self):
        dr =sqrt((self.mdat.RA-self.clusterRA)**2+(self.mdat.DEC-self.clusterDec)**2)
        galvel=self.mdat.ZDIST*300000
        vel= galvel-self.clustervel
        xlabel('dr')
        ylabel('Delta velocity')
        plot(dr,vel,'bo')


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

def allv():
    f=0
    while f<9:
        subplot(3,3,f+1)
        clusters[f].plotvel()
        f=f+1
        
