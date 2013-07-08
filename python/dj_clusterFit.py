from pylab import *
import atpy
from scipy import *
from scipy import optimize 
from LCScommon import *
import numpy as N
import binitbins as binning



class cluster:
    def __init__(self,clustername):
        infile='/home/share/research/LocalClusters/NSAmastertables/CharyElbazTables/'+clustername+'_ce_lir.fits'
        self.mdat=atpy.Table(infile)

    def plotbins(self):
        xbin,ybin,ybinerr=binning.binitbins(log10(self.mdat.LIR_ZCLUST[self.mdat.MATCHFLAG24]),(1/self.mdat.MIPS_WEIGHT[self.mdat.MATCHFLAG24]),10,min(log10(self.mdat.LIR_ZCLUST[self.mdat.MATCHFLAG24])),max(log10(self.mdat.LIR_ZCLUST[self.mdat.MATCHFLAG24])))
        figure()
        plot(xbin,ybin,'bo')
        errorbar(xbin,ybin,ybinerr,fmt=None,ecolor='b')
        xbin,ybin,ybinerr=binning.binitbins(log10(self.mdat.LIR_ZCLUST[self.mdat.MATCHFLAG24]),ones(len(self.mdat.MIPS_WEIGHT[self.mdat.MATCHFLAG24]),'f'),10,min(log10(self.mdat.LIR_ZCLUST[self.mdat.MATCHFLAG24])),max(log10(self.mdat.LIR_ZCLUST[self.mdat.MATCHFLAG24])))
        plot(xbin,ybin,'rs')
        errorbar(xbin,ybin,ybinerr,fmt=None,ecolor='r')
        gca().set_yscale('log')
        # hist(self.mdat.LIR_ZCLUST,(1/self.mdat.MIPS_WEIGHT))
        # hist(self.mdat.LIR_ZCLUST,(1/self.mdat.MIPS_WEIGHT),binthing)
           
    


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





