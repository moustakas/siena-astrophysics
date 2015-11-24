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
#import location

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
    parser.add_argument("infile2", help="Cmass data or mock data")
    parser.add_argument("--outfilename", default='default.dat', help="Outfile name")
    parser.add_argument("--range1", default=None, type=str, help="Range for first infile, input as n-n")
    parser.add_argument("--range2", default=None, type=str, help="Range for first infile, input as n-n")
    parser.add_argument('--no-plots', dest='no_plots', default=False,action='store_true', help='do not generate plots')
    parser.add_argument('--lado', dest='lado',default=False,action='store_true',help='Use Lado\'s calculations')
    parser.add_argument('--pysurvey', dest='pysurvey',default=False,action='store_true',help='Use pysurvey\'s calculations')
    parser.add_argument('--1d', dest='oned',default=False,action='store_true',help='One dimensional function')
    args=parser.parse_args()

    if args.no_plots:
        plt.switch_backend('Agg')

    infilename0 = args.infile1
    infilename1 = args.infile2

    outfilename = args.outfilename

    range1lo = None
    range1hi = None
    if args.range1 is not None:
        range1lo = int(args.range1.split('-')[0])
        range1hi = int(args.range1.split('-')[1])

    range2lo = None
    range2hi = None
    if args.range2 is not None:
        range2lo = int(args.range2.split('-')[0])
        range2hi = int(args.range2.split('-')[1])

    do_diagonal = False
    if range2lo>range1lo:
        do_diagonal = True

    # Check to see if we are using the same file for both (DD or RR)
    # or if they are different (DR)
    samefile = False
    if (infilename0==infilename1):
        samefile = True

    # Randomizing a Sample of SDSS Data
    ngals_for_calculation = 0
    nrands=0
    np.random.seed(1)

    #coords0 = jem.get_coordinates(infilename0,True,ngals_for_calculation,args.pysurvey)
    #coords1 = jem.get_coordinates(infilename1,True,nrands,args.pysurvey)

    # This is for files that do not have XYZ pre-calculated.
    coords0 = jem.get_coordinates(infilename0,False,ngals_for_calculation,args.pysurvey)
    coords1 = jem.get_coordinates(infilename1,False,nrands,args.pysurvey)
    print len(coords1)
    print len(coords0)
    print 'Read in data files and coverted to cartesian!'


    ngals0 = len(coords0)
    ngals1 = len(coords1)

    coords0cut = None
    if range1lo is not None and range1hi is not None:
        coords0cut = coords0[range1lo:range1hi]
    else:
        coords0cut = coords0

    coords1cut = None
    if range2lo is not None and range2hi is not None:
        coords1cut = coords1[range2lo:range2hi]
    else:
        coords1cut = coords1
   
    # This is for the histogram.
    nbins=200
    rangeval=200

    if args.oned:
        nbins*=2

    tot_freq = np.zeros((nbins,nbins)) 
    if args.oned==True:
        tot_freq = np.zeros(nbins) 

    
    # Figure out the chunking.
    

    chunk_size = 50
    nchunks = len(coords0cut)/chunk_size     #ngals_for_calculation/chunk_size

    ncalcs_per_chunk = chunk_size*len(coords1cut) #chunk_size*ngals1

    # These will store the calculations.
    paras = np.zeros(ncalcs_per_chunk)
    perps = np.zeros(ncalcs_per_chunk)

    indexlo = 0
    indexhi = 0

    #Calculation Loop
    for j in xrange(nchunks):
        lo = j*chunk_size
        hi = (j+1)*chunk_size
        #print "Performing calculations for DR %d chunk: %d-%d" % (j,lo,hi)

        paras *= 0.
        perps *= 0.

        #for i,r0 in enumerate(coords0[lo:hi]):
        for i in range(lo,hi):
            r0 = coords0cut[i]

            lo1 = 0
            if samefile:
                lo1 = i
                if do_diagonal==False:
                    lo1 += 1

            indexhi += len(coords1cut[lo1:])

            other_gals = coords1cut[lo1:]

            # Calc 1D using the pysurvey distance calc.
            if args.lado==False and args.oned==False and args.pysurvey==True:
                if len(other_gals)>0:
                    temp_paras,temp_perps = jem.pysurvey_distance(r0,other_gals)

            # Calc para and perp ``our" way. 
            elif args.lado==False and args.oned==False and args.pysurvey==False:
                temp_paras,temp_perps = jem.our_para_perp(r0,other_gals)

            # Calc Lado's way
            elif args.lado==True and args.oned==False and args.pysurvey==False:
                temp_paras,temp_perps = jem.lado_para_perp(r0,other_gals)

            # Calc just the 1D
            elif args.lado==False and args.oned==True and args.pysurvey==False:
                temp_paras,temp_perps = jem.one_dimension(r0,other_gals)

            paras[indexlo:indexhi] = temp_paras
            perps[indexlo:indexhi] = temp_perps

            indexlo = indexhi
        
        # Histogram the values.
        if args.oned==False:
            hist=plt.hist2d(perps[0:indexhi],paras[0:indexhi],bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
        else:
            # THE NUMPY HISTOGRAM SEEMS TO BE FASTER THAN THE MATPLOTLIB ONE HERE. 
            hist=np.histogram(paras[0:indexhi],bins=nbins,range=(0,rangeval))

        # Mirror the negative perps
        #hist=plt.hist2d(-1*perps[0:indexhi],paras[0:indexhi],bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
        tot_freq += hist[0]

        indexlo=0
        indexhi=0

        del hist

        print tot_freq.sum()
   
    print 'Point:'
    #print tot_freq[199,101]

    # Do this for the 1D
    xvals = np.linspace(0,rangeval,nbins+1) # THIS IS PROBABLY NOT EXACTLY CORRECT

    if args.no_plots==False:
        print 'Final Plot'    
        if args.oned==False:
            extent = [-rangeval,rangeval, -rangeval,rangeval]
            fig = plt.figure()
            axes = fig.add_subplot(1,1,1)
            print 'Imshow'
            print tot_freq
            ret = axes.imshow(tot_freq,extent=extent,interpolation='nearest') #,origin=origin,cmap=cmap,axes=axes,aspect=aspect
        else:
            plt.plot(xvals[0:nbins],tot_freq,'k.')

        plt.show()

    if args.oned==False:
        print('Writing {}'.format(outfilename))
        np.savetxt(outfilename,tot_freq)
    else:
        print('Writing {}'.format(outfilename))
        outfile = open(outfilename,"w")
        output = "%d,%d,0,0\n" % (ngals0,ngals1)
        outfile.write(output)
        for i in xrange(nbins):
            output = "%f,%f,%f,%f\n" % (xvals[i],(xvals[i+1]+xvals[i])/2.,xvals[i+1],tot_freq[i])
            outfile.write(output)
        outfile.close()




################################################################################
if __name__=='__main__':
    main()

