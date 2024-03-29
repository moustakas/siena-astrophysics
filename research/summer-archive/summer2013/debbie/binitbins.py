import pylab as pl
import atpy
import numpy as np

def binitbins(x,w,nbin,xmin,xmax):#use equally spaced bins
    if len(x) != len(w):
        print 'array and weights must be same length'
        return
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
        try:
            ybin[i]=sum(xscaled)
            ybinerr[i]=np.sqrt(ybin[i])
        except ZeroDivisionError:
            ybin[i]=0
            ybinerr[i]=0
            
    return xbin,ybin,ybinerr
