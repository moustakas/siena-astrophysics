from pylab import *
import atpy
from scipy import *
from scipy import optimize 
from LCScommon import *
import numpy as N


class cluster:
    def __init__(self,clustername):
        infile='/home/share/research/LocalClusters/NSAmastertables/'+clustername+'_NSAmastertable.fits'
        self.mdat=atpy.Table(infile)

    def binitbins(xmin,xmax,nbin,x,y):#use equally spaced bins
        dx=float((xmax-xmin)/(nbin))
        xbin=N.arange(xmin,(xmax),dx)+dx/2.
        ybin=N.zeros(len(xbin),'d')
        ybinerr=N.zeros(len(xbin),'d')
        xbinnumb=N.array(len(self.mdat.ABSMAG[:3]),'d')
        x1=N.compress((x >= xmin) & (x <= xmax),x)
        y1=N.compress((x >= xmin) & (x <= xmax),y) 
        x=x1
        y=y1
        xbinnumb=((x-xmin)*nbin/(xmax-xmin))#calculate x  bin number for each point
        j=-1
        for i in range(len(xbin)):
            ydata=N.compress(abs(xbinnumb-float(i))<.5,y)
            try:
                ybin[i]=pylab.median(ydata)
                ybinerr[i]=pylab.std(ydata)/N.sqrt(float(len(ydata)))
            except ZeroDivisionError:
                ybin[i]=0.
                ybinerr[i]=0.
                return xbin,ybin,ybinerr
    

    def lum(self):
        t=hist(self.mdat.ABSMAG[:,3],bins=100)
        clf()
        bins=t[1]
        yvalue=t[0]
        yerror=sqrt(t[0])
        bincenter=[]
        for x in range(len(yvalue)):
            bincenter.append((bins[x]+(bins[x+1]))/2)
        yplot=log10(yvalue)
        plot(bincenter,yplot,'ko')
        yerrup=log10(yvalue+yerror)-log10(yvalue)
        yerrdown=log10(yvalue)-log10(yvalue-yerror)
        yerrboth=zip(yerrdown,yerrup)
        yerrboth=transpose(yerrboth)
        errorbar(bincenter,yplot,yerrboth)

        logphistar=0
        alpha=-1
        loglstar=.3

       
        fitfunc = lambda  x,: bincenter[0]*yvalue+ bincenter[1] * x
        errfunc = lambda  x,yerror,logphistar,alpha,loglstar: (((log(10)*10.**logphistar)*((10.**array(bincenter))/(10.**loglstar))**(alpha+1.)*(exp(-(10.**array(bincenter))/(10.**loglstar)))) - fitfunc(x)) / yerror
        pinit = [1.0, -1.0]
        out = optimize.leastsq(errfunc(x,yerror,logphistar,alpha,loglstar), pinit,args=(bincenter))
        print errfunc    

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
        xlabel('g band log(L/L$_\odot$)')
        ylabel('Number of Galaxies')
        (clusters[f-1]).lum()
        legend([clusternames[f-1]])
        Lvalues=linspace(-2.,1.,100)
        Llin=10**(Lvalues)
        f=f+1
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
    pinit = [1.0, -1.0]
    out = optimize.leastsq(errfunc, pinit,args=(logx, logy, logyerr), full_output=1)



def binitbins(xmin,xmax,nbin,x,y):#use equally spaced bins
    dx=float((xmax-xmin)/(nbin))
    xbin=N.arange(xmin,(xmax),dx)+dx/2.
    ybin=N.zeros(len(xbin),'d')
    ybinerr=N.zeros(len(xbin),'d')
    xbinnumb=N.array(len(x),'d')
    x1=N.compress((x >= xmin) & (x <= xmax),x)
    y1=N.compress((x >= xmin) & (x <= xmax),y) 
    x=x1
    y=y1
    xbinnumb=((x-xmin)*nbin/(xmax-xmin))#calculate x  bin number for each point 
    j=-1
    for i in range(len(xbin)):
	ydata=N.compress(abs(xbinnumb-float(i))<.5,y)
	try:
	    #ybin[i]=N.average(ydata)
	    ybin[i]=pylab.median(ydata)
	    ybinerr[i]=pylab.std(ydata)/N.sqrt(float(len(ydata)))
	except ZeroDivisionError:
	    ybin[i]=0.
	    ybinerr[i]=0.
    return xbin,ybin,ybinerr
