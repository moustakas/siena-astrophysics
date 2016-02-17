import astropy.io 
from astropy.io import fits
import numpy as np
import sys
import random
import math
import scipy
import scipy.spatial
import matplotlib.pylab as plt
from astropy.cosmology import FlatLambdaCDM
import time
import argparse
import timeit

import jem_utilities as jem

def main():
    """Given command-line arguments, will return
        frequency arrays for galactic distances.

    Args:
        See Below

    Returns:
        Dependant on command-line arguments.
    """
    parser= argparse.ArgumentParser()
    parser.add_argument("infile1", help="Cmass data or mock data")
    parser.add_argument("infile2", default=None, help="Cmass data or mock data")
    parser.add_argument("--outfilename",default=None, help="Outfile name to save data")
    parser.add_argument("--range1", default=None, type=str, help="Range for first infile, input as n-n")
    parser.add_argument("--range2", default=None, type=str, help="Range for first infile, input as n-n")
    parser.add_argument('--pysurvey', dest='pysurvey',default=False,action='store_true',help='Use pysurvey\'s calculations')
    parser.add_argument('--1d', dest='oned',default=False,action='store_true',help='One dimensional function')
    parser.add_argument('--maxdist',type=int,default=200,help='Maximum Distance,Default=200Mpc')
    parser.add_argument('--distres',type=int,default=1,help='Resolution,Default=1Mpc/bin')
    args=parser.parse_args()

  
    
    start_time=timeit.timeit()

    infilename1=args.infile1
    infilename2=args.infile2

    if infilename2 is None:
        infilename2=infilename1
    
    range1=args.range1
    range2=args.range2

    rangeval=args.maxdist
    distres=args.distres
    nbins=rangeval/distres
    print rangeval                    
    print distres
    print nbins

    DD_calc=jem.twopoint_hist(infilename1,infilename1,
                         nbins,rangeval,range1,range2,
                                       oned=args.oned,)

    
    outfilenameDD =  'DD.dat'
    np.savetxt(outfilenameDD,DD_calc)
    
    RR_calc=jem.twopoint_hist(infilename2,infilename2,
                         nbins,rangeval,range1,range2,
                                      oned=args.oned)
                                      

    
    outfilenameRR =  'RR.dat'
    np.savetxt(outfilenameRR,RR_calc)
    
    DR_calc=jem.twopoint_hist(infilename1,infilename2,
                         nbins,rangeval,range1,range2,
                                       oned=args.oned)
                                       
        
        
    outfilenameDR = 'DR.dat'    
    np.savetxt(outfilenameDR,DR_calc)
    
    Xi=jem.corr_est(DD_calc,DR_calc,RR_calc,2000,2000,nbins,oned=args.oned)
    

    jem.corr_plot(Xi,-rangeval,rangeval,-rangeval,rangeval,"Title","Xlabel","Ylabel",oned=args.oned)

    #Saving
    outfilename=args.outfilename
    if outfilename:
        print 'Saving'
        outfilenameDD = outfilename + 'DD.dat'
        outfilenameDR = outfilename + 'DR.dat'
        outfilenameRR = outfilename + 'RR.dat'
        np.savetxt(outfilenameDD,DD_calc)
        np.savetxt(outfilenameDR,DR_calc)
        np.savetxt(outfilenameRR,RR_calc)
    
    end_time=timeit.timeit()
    sec=end_time-start_time
    hour=sec/360
    remmin=sec%hour
    remsec=sec%remmin

    print "This code took %f hours, %f minutes, and %f seconds to run" % (hour, remmin, remsec)

if __name__ == "__main__":
    main()
