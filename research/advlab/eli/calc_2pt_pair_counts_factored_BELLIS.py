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

    # Check to see if we are using the same file for both (DD or RR)
    # or if they are different (DR)
    samefile = False
    if (infilename0==infilename1):
        samefile = True

    # Assume reading in from a text file.
    # Return ra,dec (in radians), redshift z, and weight.
    start = time.time()
    print "Reading in data...."
    print "Opening ",infilename0
    coords0 = jem.get_coordinates_with_weight(infilename0)
    print "Opening ",infilename1
    coords1 = jem.get_coordinates_with_weight(infilename1)
    print "Read in data......"
    time_read = time.time()
    print "Time to read in data %f" % (time_read - start)
    print "Total execution time %f" % (time.time() - start)

    # Convert to x,y,z
    print "Converting to x,y,z...."
    coords0 = jem.radecredshift2xyz_with_weights(coords0)
    coords1 = jem.radecredshift2xyz_with_weights(coords1)
    time_convert = time.time()
    print "Time to convert data %f" % (time_convert - time_read)
    print "Total execution time %f" % (time.time() - start)

    ngals0 = len(coords0)
    ngals1 = len(coords1)

    print "Ngals 0/1: ",ngals0,ngals0

    ############################################################################
    # Break up into voxels.
    ############################################################################

    maxsep=200
    nbins=20
    #'''
    print "Breaking up into voxels..."
    # ngrids tells you how many galaxies are along each axis. 
    voxels0,voxels1,ngrids,gridwidths,loranges,hiranges  = jem.voxelize_the_data(coords0,coords1,maxsep=maxsep)
    time_vox = time.time()
    print "Time to voxelize data %f" % (time_vox - time_convert)
    print "Total execution time %f" % (time.time() - start)

    print "Performing the pair counts...."
    pair_counts = jem.do_pair_counts(voxels0,voxels1,ngrids,nbins=nbins,maxrange=maxsep,samefile=samefile)
    #pair_counts = jem.do_pair_counts_2d(voxels0,voxels1,ngrids,nbins=nbins,maxrange=maxsep,samefile=samefile)
    time_pc = time.time()
    print "PAIR COUNTS"
    print  pair_counts.shape
    print "Time to perform pair counts %f" % (time_pc - time_pc)
    print "Total execution time %f" % (time.time() - start)
    #'''

    '''
    # Alternatively (for small numbers, <10k)
    coords0t=np.column_stack((coords0[:,0],coords0[:,1],coords0[:,2]))
    coords1t=np.column_stack((coords1[:,0],coords1[:,1],coords1[:,2]))
    print len(coords0t)
    print len(coords0t[0])
    #distances=scipy.spatial.distance.pdist(coords)
    distances=scipy.spatial.distance.cdist(coords0t,coords1t)
    #print np.sort(distances[distances<maxsep])
    hist=np.histogram(distances,bins=nbins,range=(0,maxsep))
    pair_counts = hist[0]
    '''

    # New way to normalize the weighting
    tot_weight2 = 0.
    w0 = coords0.transpose()[3]
    w1 = coords1.transpose()[3]
    for i in w0:
       tot_weight2 += (i*w1).sum()

    #'''
    if samefile==True:
        tot_weight2 -= (w0*w1).sum()
        tot_weight2 /= 2.
        tot_weight3 = (w0*w1).sum()
    else:
        tot_weight3 = 1.0 # Otherwise, not necessarily the same size.
    #'''

    tot_weight0 = w0.sum()
    tot_weight1 = w1.sum()

    print "Tot weight calc the new way: %f %f %f %f" % (tot_weight2,tot_weight0,tot_weight1,tot_weight3)

    #pair_counts /= tot_weight


    #print "Sum: ",sum(pair_counts)


    # Do this for the 1D
    xvals = np.linspace(0,maxsep,nbins+1) # THIS IS PROBABLY NOT EXACTLY CORRECT

    if args.no_plots==False:
        print 'Final Plot'    
        plt.plot(xvals[0:nbins],pair_counts,'k.')

        plt.show()

    if args.oned==False:
        print('Writing {}'.format(outfilename))
        outfile = open(outfilename,"w")
        line_one = np.zeros(nbins)
        line_one = np.insert(line_one,0,ngals0)
        line_one = np.insert(line_one,1,ngals1)
        output = line_one
        outfile.write(output)
        output = pair_counts
        print output
        outfile.write(output)
        outfile.close()

    else:
        print('Writing {}'.format(outfilename))
        outfile = open(outfilename,"w")
        output = "%d,%d,%f,%f,%f,%f\n" % (ngals0,ngals1,tot_weight0,tot_weight1,tot_weight2,tot_weight3)
        outfile.write(output)
        for i in xrange(nbins):
            output = "%f,%f,%f,%f,0,0\n" % (xvals[i],(xvals[i+1]+xvals[i])/2.,xvals[i+1],pair_counts[i])
            print output.rstrip()
            outfile.write(output)
        outfile.close()

    print "Finished writing out data."
    print "Total execution time %f" % (time.time() - start)


################################################################################
if __name__=='__main__':
    main()

