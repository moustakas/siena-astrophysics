import pylab as pl
import atpy
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
<<<<<<< HEAD
        
=======
>>>>>>> ac8aa156f00f83c60b298d488f2aff6313a9fac3
