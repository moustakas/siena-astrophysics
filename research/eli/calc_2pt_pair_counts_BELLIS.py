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
    #coords0 = jem.get_coordinates(infilename0,False,ngals_for_calculation,args.pysurvey)
    #coords1 = jem.get_coordinates(infilename1,False,nrands,args.pysurvey)
    #print len(coords1)
    #print len(coords0)
    #print 'Read in data files and coverted to cartesian!'

    # Trying something based on conversation with Rose. 
    # Return ra,dec,z in radians.
    coords0 = jem.get_coordinates(infilename0,False,0,return_radecz=True)
    coords1 = jem.get_coordinates(infilename1,False,0,return_radecz=True)

    # Convert to x,y,z
    print "Converting to x,y,z...."
    coords0 = jem.radecredshift2xyz(coords0[:,0],coords0[:,1],coords0[:,2])
    coords1 = jem.radecredshift2xyz(coords1[:,0],coords1[:,1],coords1[:,2])

    print coords0
    #exit()

    ngals0 = len(coords0)
    ngals1 = len(coords1)

    print ngals0
    print ngals0
    #print 'Read in data files and left as ra,dec, and redshift!'


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
   
    ############################################################################
    # Break up into voxels.
    ############################################################################
    loranges = [min(coords0cut[:,0]),min(coords0cut[:,1]),min(coords0cut[:,2])]
    hiranges = [max(coords0cut[:,0]),max(coords0cut[:,1]),max(coords0cut[:,2])]
    ngrids,gridwidths = jem.define_ranges(loranges,hiranges, maxsep=200)

    grid_coords0 = jem.assign_grid_coordinate(coords0cut, loranges, hiranges, gridwidths)
    grid_coords1 = jem.assign_grid_coordinate(coords1cut, loranges, hiranges, gridwidths)

    print ngrids
    print gridwidths
    print grid_coords0
    #exit()

    # Subdivide into voxels.

    # First make our holder of coordinates, broken up by voxel.
    # Initialize them
    voxels0 = []
    for i in range(0,ngrids[0]):
        voxels0.append([])
        for j in range(0,ngrids[1]):
            voxels0[i].append([])
            for k in range(0,ngrids[2]):
                voxels0[i][j].append([])

    voxels1 = []
    for i in range(0,ngrids[0]):
        voxels1.append([])
        for j in range(0,ngrids[1]):
            voxels1[i].append([])
            for k in range(0,ngrids[2]):
                voxels1[i][j].append([])

    # Fill them
    for i in range(0,ngals0):
        a,b,c = grid_coords0[0][i],grid_coords0[1][i],grid_coords0[2][i]
        #print a,b,c
        #print a,b,c,coords0[i],loranges,hiranges
        voxels0[a][b][c].append(coords0[i])
        #print voxels0
    #print voxels0[4][10][5]
    #exit()

    for i in range(0,ngals1):
        a,b,c = grid_coords1[0][i],grid_coords1[1][i],grid_coords1[2][i]
        voxels1[a][b][c].append(coords1[i])

    #print voxels1
    # Convert lists to arrays
    for i in range(0,ngrids[0]):
        for j in range(0,ngrids[1]):
            for k in range(0,ngrids[2]):
                voxels0[i][j][k] = np.array(voxels0[i][j][k])

    for i in range(0,ngrids[0]):
        for j in range(0,ngrids[1]):
            for k in range(0,ngrids[2]):
                voxels1[i][j][k] = np.array(voxels1[i][j][k])

    #exit()



    ############################################################################
    # This is for the histogram.
    nbins=10
    rangeval=200

    if args.oned:
        nbins*=2

    tot_freq = np.zeros((nbins,nbins)) 
    if args.oned==True:
        tot_freq = np.zeros(nbins) 

    
    #Calculation Loop
    for ii in range(0,ngrids[0]):
        for jj in range(0,ngrids[1]):
            for kk in range(0,ngrids[2]):

                c0 = voxels0[ii][jj][kk] 

                iimax = ii+2
                if iimax>=ngrids[0]:
                    iimax = ngrids[0]
                jjmax = jj+2
                if jjmax>=ngrids[1]:
                    jjmax = ngrids[1]
                kkmax = kk+2
                if kkmax>=ngrids[2]:
                    kkmax = ngrids[2]

                for aa in range(ii,iimax):
                    for bb in range(jj,jjmax):
                        for cc in range(kk,kkmax):

                            c1 = voxels1[aa][bb][cc] 

                            #print "NVALS IN VOXELS: ",len(c0),len(c1)

                            if len(c0)>0 and len(c1)>0:
                                print ii,jj,kk, aa, bb, cc, len(c0),len(c1)

                            if len(c0)==0 or len(c1)==0:
                                continue

                            # These will store the calculations.
                            for index,r0 in enumerate(c0):
                                '''
                                for i in range(lo,hi):
                                    r0 = c0[i]

                                    lo1 = 0
                                    if samefile:
                                        lo1 = i
                                        if do_diagonal==False:
                                            lo1 += 1
                                    '''

                                if samefile and aa==ii and bb==jj and cc==kk:
                                    paras,perps = jem.one_dimension(r0,c1[index+1:])
                                else:
                                    paras,perps = jem.one_dimension(r0,c1)

                                hist=np.histogram(paras,bins=nbins,range=(0,rangeval))

                                tot_freq += hist[0]

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

