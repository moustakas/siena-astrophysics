import pylab as pl
import atpy
import scipy 
from scipy import optimize 
from LCScommon import *
import numpy as np

def binitbins(x,w,nbin,xmin,xmax):#use equally spaced bins
    #if self.mdat.MATCHFLAG24 == True:
    dx=float((xmax-xmin)/(nbin))
    
    xbin=np.arange(xmin,(xmax),dx)+dx/2.
    
    ybin=np.zeros(len(xbin),'d')
    print ybin
    ybinerr=np.zeros(len(xbin),'d')
    print ybinerr
    xbinnumb=np.array(len(x),'d')
    
    x=x[((x >= xmin) & (x <= xmax))]

    
    xbinnumb=((x-xmin)*nbin/(xmax-xmin))-1#calculate x  bin number for each point
    print x,xbinnumb,xbin
    for i in range(len(xbin)):
        xbindata=x[(abs(xbinnumb-float(i))<.5)]
        wbindata=w[(abs(xbinnumb-float(i))<.5)]
        xscaled=xbindata*wbindata
             
                   
        #print abs(xbinnumb-float(i))
        #print "ydata"
        try:
            ybin[i]=sum(xscaled)
            ybinerr[i]=sqrt(xbindata)*wbindata
    return xbin,ybin,ybinerr
        


class cluster:
    def __init__(self,clustername):
        infile='/home/share/research/LocalClusters/NSAmastertables/CharyElbazTables/'+clustername+'_ce_lir.fits'
        self.mdat=atpy.Table(infile)
        





mkw11=cluster('MKW11')
mkw8=cluster('MKW8')
coma=cluster('Coma')
herc=cluster('Hercules')
ngc=cluster('NGC6107')
awm4=cluster('AWM4')
a1367=cluster('A1367')
a2052=cluster('A2052')
a2063=cluster('A2063')
