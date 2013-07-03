#!/usr/bin/env python
from pylab import *
import atpy
from scipy import optimize



class ediscs:
    def __init__(self):
        infile='/home/share/research/ediscs/ediscs_ircatalogs_spec.fits'
        self.ircats=atpy.Table(infile)
        flag1=(self.ircats.MATCHFLAG24==1)
        flag2=(self.ircats.SPECMEMBFLAG==1)
        mflag=where(flag1)
        sflag=where(flag2)
        superflag=flag1 & flag2
        #edimatch=self.ircats[mflag]
        #edispec=self.ircats[sflag]
        detections=self.ircats[superflag]
        print "Print Detections"
        print detections
    
    def ediscshist(self):
        hist(detections, bins=100)
        figure()
    #def __init__(self):
        #infile='/home/

    

    
edi=ediscs()
