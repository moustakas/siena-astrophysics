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
    print "Reading in data...."
    coords0 = jem.get_coordinates_with_weight(infilename0)
    coords1 = jem.get_coordinates_with_weight(infilename1)
    print "Read in data......"

    # Convert to x,y,z
    print "Converting to x,y,z...."
    coords0 = jem.radecredshift2xyz_with_weights(coords0)
    coords1 = jem.radecredshift2xyz_with_weights(coords1)
    print "Converted to x,y,z...."

    ngals0 = len(coords0)
    ngals1 = len(coords1)

    print ngals0
    print ngals0
    #print 'Read in data files and left as ra,dec, and redshift!'


    ############################################################################
    # Break up into voxels.
    ############################################################################

    maxsep=200
    nbins=20
    # ngrids tells you how many galaxies are along each axis. 
    #'''
    voxels0,voxels1,ngrids,gridwidths,loranges,hiranges  = jem.voxelize_the_data(coords0,coords1,maxsep=maxsep)

    # Count the entries as a sanity check.
    tot = 0
    for i in voxels0:
        for j in i:
            for k in j:
                tot += len(k)

    print "Tot: ",tot
    print ngrids
    #exit()

    pair_counts = jem.do_pair_counts(voxels0,voxels1,ngrids,nbins=nbins,maxrange=maxsep,samefile=samefile)

    print pair_counts
    #exit()
    #'''


    '''
    # Alternatively (for small numbers, <10k)
    coords=np.column_stack((coords0[:,0],coords0[:,1],coords0[:,2]))
    print len(coords)
    print len(coords[0])
    distances=scipy.spatial.distance.pdist(coords)
    #print np.sort(distances[distances<maxsep])
    hist=np.histogram(distances,bins=nbins,range=(0,maxsep))
    pair_counts = hist[0]
    '''


    print "Sum: ",sum(pair_counts)


    # Do this for the 1D
    xvals = np.linspace(0,maxsep,nbins+1) # THIS IS PROBABLY NOT EXACTLY CORRECT

    if args.no_plots==False:
        print 'Final Plot'    
        plt.plot(xvals[0:nbins],pair_counts,'k.')

        plt.show()

    print('Writing {}'.format(outfilename))
    outfile = open(outfilename,"w")
    output = "%d,%d,0,0\n" % (ngals0,ngals1)
    outfile.write(output)
    for i in xrange(nbins):
        output = "%f,%f,%f,%f\n" % (xvals[i],(xvals[i+1]+xvals[i])/2.,xvals[i+1],pair_counts[i])
        print output.rstrip()
        outfile.write(output)
    outfile.close()




################################################################################
if __name__=='__main__':
    main()

