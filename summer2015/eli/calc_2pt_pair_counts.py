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

################################################################################
# Get magnitude of a vector
################################################################################
def mag(vec):

    m = None
    if type(vec)==np.ndarray:
        m = np.sqrt(vec[:,0]**2 + vec[:,1]**2 + vec[:,2]**2)
    else:
        m = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

    return m

################################################################################

################################################################################
# Converting RA and Dec and redshift to Cartesian
################################################################################
def radecredshift2xyz(ra,dec,redshift):

    # Comoving Distances In Mpc
    cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
    comdist=cosmo.comoving_distance(redshift).value

    # Convert spherical to Cartesian Coordinates
    x=comdist*np.sin(dec)*np.cos(ra)
    y=comdist*np.sin(dec)*np.sin(ra)
    z=comdist*np.cos(dec)

    coords=np.column_stack((x,y,z))

    return coords

################################################################################
################################################################################


#@profile
'''
def pair_counts():
    #t0=time.time()
'''

################################################################################
def get_coordinates(infilename,maxgals=0):

    isdatafile = False
    if(infilename.find('fits')>=0):
        isdatafile = True

    if isdatafile:
        # Opening FITS file (SDSS Data) 'a' Argument = dr10cmassnorth.fits
        print 'Reading in FITS Data'
        hdulist1=fits.open(infilename)
        hdulist1.info()
        h1=hdulist1[1]
        data=h1.data
        del h1


        # Radians
        ra=(data['PLUG_RA'])*((math.pi)/180)
        dec=((math.pi)/2)-((data['PLUG_DEC'])*((math.pi)/180))
        redshift = data['Z']

        del data
    else:
        # Opening txt file (Mocks) 'b'
        print 'Reading in Text File'
        r=np.loadtxt(infilename)

        # Radians
        ra=(r[:,0])*((math.pi)/180)
        dec=((math.pi)/2)-((r[:,1])*((math.pi)/180))
        redshift=r[:,2]

    # Made some common cuts
    index0 = redshift<0.7
    index1 = redshift>0.43
    index = index0*index1

    ra = ra[index]
    dec = dec[index]
    redshift = redshift[index]

    # Grab a subsample if we only want a few galaxies
    if maxgals>0:

        a=np.arange(0,len(redshift))
        np.random.shuffle(a)

        ra=ra[a[0:maxgals]]
        dec=dec[a[0:maxgals]]
        redshift=redshift[a[0:maxgals]]
        del a

    coords = radecredshift2xyz(ra,dec,redshift)

    del ra,dec,redshift

    return coords

################################################################################



################################################################################
def main():

    parser= argparse.ArgumentParser()
    parser.add_argument("infile1", help="Cmass data or mock data")
    parser.add_argument("infile2", help="Cmass data or mock data")
    parser.add_argument("--outfilename", default='default.dat', help="Outfile name")
    parser.add_argument("--range1", default=None, type=str, help="Range for first infile, input as n-n")
    parser.add_argument("--range2", default=None, type=str, help="Range for first infile, input as n-n")
    parser.add_argument('--no-plots', dest='no_plots', default=False,action='store_true', help='do not generate plots')
    parser.add_argument('--lado', dest='lado',default=False,action='store_true',help='Use Lado\'s calculations')
    #parser.add_argument('--oned', dest='1d' ,default=False,action='store_true',help='One dimensional Function')
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

    coords0 = get_coordinates(infilename0,ngals_for_calculation)
    coords1 = get_coordinates(infilename1,nrands)
    
    print 'Read in data files and coverted to cartesian!'

    ################################################################################


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
   
    chunk_size = 50
    nchunks = len(coords0cut)/chunk_size     #ngals_for_calculation/chunk_size
    nbins=200
    rangeval=300

    tot_freq = np.zeros((nbins,nbins)) 

    ncalcs_per_chunk = chunk_size*len(coords1cut) #chunk_size*ngals1

    paras = np.zeros(ncalcs_per_chunk)
    perps = np.zeros(ncalcs_per_chunk)

    indexlo = 0
    indexhi = 0
    ################################ ONE DIMENSION ##########################################

    #if args.1d:
     #   x0=coords0cut[:,0]
      ## y0=coords0cut[:,1]
        #y1=coords1cut[:,1]
        #z0=coords0cut[:,2]
        #z1=coords1cut[:,2]
        #lo1=0
        #if samefile:
            
        #for j in xrange(nchunks):
         #   lo = j*chunk_size
          #  hi = (j+1)*chunk_size
           # for i in range(lo,hi):
            #    val=[np.sqrt(((x0[i]-x1[i+1:]**2)+(y1[i]-
                
                
                
            
        
            
    ################################## LADO'S CODE ############################
    if args.lado:
        
        for j in xrange(nchunks):
            lo = j*chunk_size
            hi = (j+1)*chunk_size
            #print "Performing calculations for DR %d chunk: %d-%d" % (j,lo,hi)
            print j
            paras *= 0.
            perps *= 0.

            #for i,r0 in enumerate(coords0[lo:hi]):
            for i in range(lo,hi):
                x0=coords0cut[i,0]
                y0=coords0cut[i,1]
                z0=coords0cut[i,2]
                lo1 = 0
                if samefile:
                    lo1 = i
                    if do_diagonal==False:
                        lo1 += 1

                indexhi += len(coords1cut[lo1:])
        
                
                #x1=coords1cut[lo1:,0]
                #y1=(coords1cut[lo1:,1])
                #z1=(coords1cut[lo1:,2])

                #x1=coords1cut[lo1:,0]
                #y1=coords1cut[lo1:,1]
                #z1=coords1cut[lo1:,2]
                

                d1 = x0-(coords1cut[lo1:,0])
                #print len(d1)#(coords1cut[lo1:,0])
                d2 = y0-(coords1cut[lo1:,1]) #(coords1cut[lo1:,1])
                d3 = z0-(coords1cut[lo1:,2]) #(coords1cut[lo1:,2])
                
                r2l = d1*d1 + d2*d2 + d3*d3
                gd1 = np.sqrt((x0**2)+(y0**2)+(z0**2))
                gd2 = np.sqrt(((coords1cut[lo1:,0])**2)+((coords1cut[lo1:,1])**2)+((coords1cut[lo1:,2])**2))
                rat = gd1/gd2

                xb = x0 + (coords1cut[lo1:,0])*rat
                yb = y0 + (coords1cut[lo1:,1])*rat
                zb = z0 + (coords1cut[lo1:,2])*rat

                db2 = xb*xb + yb*yb + zb*zb

                mu = np.absolute(((xb*d1 + yb*d2 + zb*d3)/np.sqrt(r2l)/np.sqrt(db2)))
                rr = np.sqrt(r2l)
                 
                rpar=rr*mu
                rperp=rr*np.sqrt(1-(mu*mu))
                #print indexlo,indexhi,type(rpar)
                
                paras[indexlo:indexhi] = rpar
                perps[indexlo:indexhi] = rperp    

                indexlo = indexhi

            hist=plt.hist2d(perps[0:indexhi],paras[0:indexhi],bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
            tot_freq += hist[0]
            
            # Mirror the negative perps
            #hist=plt.hist2d(-1*perps[0:indexhi],paras[0:indexhi],bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
            #tot_freq += hist[0]
            #hist=plt.hist2d(perps[0:indexhi],-1*paras[0:indexhi],bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
            #tot_freq += hist[0]
            #hist=plt.hist2d(-1*perps[0:indexhi],-1*paras[0:indexhi],bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
            #tot_freq += hist[0]
            #print type(hist1[0])
            #frequ1=hist1[0]
            #plt.close()

            indexlo=0
            indexhi=0
            #del paras
            #del perps
            del hist

            #print tot_freq
            #tot_freq[(nbins/2),(nbins/2)]=0
            print tot_freq.sum()

        if args.no_plots==False:
            print 'Final Plot'    
            extent = [-rangeval,rangeval, -rangeval,rangeval]
            fig = plt.figure()
            axes = fig.add_subplot(1,1,1)
            print 'Imshow'
            print tot_freq
            ret = axes.imshow(tot_freq,extent=extent,interpolation='nearest') #,origin=origin,cmap=cmap,axes=axes,aspect=aspect
            plt.show()

        print('Writing {}'.format(outfilename))
        print 'LADO'
        np.savetxt(outfilename,tot_freq)



####################################OUR CODE############################

    else:        
        for j in xrange(nchunks):
            lo = j*chunk_size
            hi = (j+1)*chunk_size
            #print "Performing calculations for DR %d chunk: %d-%d" % (j,lo,hi)

            paras *= 0.
            perps *= 0.

            #for i,r0 in enumerate(coords0[lo:hi]):
            for i in range(lo,hi):
                r0 = coords0cut[i]
                #print i,lo,hi
                #print lo+i
                #print r0

                
                lo1 = 0
                if samefile:
                    lo1 = i
                    if do_diagonal==False:
                        lo1 += 1

                indexhi += len(coords1cut[lo1:])
                

                # First compute R_LOS and dR
                R_LOS = (r0 + coords1cut[lo1:])/2
                dR = coords1cut[lo1:] - r0
                R_LOS_mag = mag(R_LOS)

                # Dot product
                
                R_para = (dR[:,0]*R_LOS[:,0] + dR[:,1]*R_LOS[:,1] + dR[:,2]*R_LOS[:,2])/R_LOS_mag
                
                dR_mag = mag(dR)
                
                # Make use of the Pythagorean theorem
                R_perp = np.sqrt(dR_mag*dR_mag - R_para*R_para)
                
                #print i,lo1,indexlo,indexhi,len(R_para),len(paras)

                paras[indexlo:indexhi] = R_para
                perps[indexlo:indexhi] = R_perp

                #print len(paras[paras!=0])

                indexlo = indexhi
                #nperps1 += negR_perp1.tolist()
                #if i%(chunk_size/4)==0:
                    #print i
            

            #print len(paras)
            #print len(perps)
            #newperps1=np.concatenate((perps1,nperps1))
            #newparas1len=np.concatenate((paras1,paras1))

            #print 'Histogram1'

            #print paras[0:10]
            #print indexhi
            hist=plt.hist2d(perps[0:indexhi],paras[0:indexhi],bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
            tot_freq += hist[0]
            
            # Mirror the negative perps
            hist=plt.hist2d(-1*perps[0:indexhi],paras[0:indexhi],bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
            tot_freq += hist[0]

            #print type(hist1[0])
            #frequ1=hist1[0]
            #plt.close()

            indexlo=0
            indexhi=0
            #del paras
            #del perps
            del hist

            #print tot_freq
            #tot_freq[(nbins/2),(nbins/2)]=0
            print tot_freq.sum()
       
        #tot_freq[(nbins/2),(nbins/2)]=0
        if args.no_plots==False:
            print 'Final Plot'    
            extent = [-rangeval,rangeval, -rangeval,rangeval]
            fig = plt.figure()
            axes = fig.add_subplot(1,1,1)
            print 'Imshow'
            ret = axes.imshow(tot_freq,extent=extent,interpolation='nearest') #,origin=origin,cmap=cmap,axes=axes,aspect=aspect
            plt.show()

        print('Writing {}'.format(outfilename))
        print 'OUR VERSION'
        np.savetxt(outfilename,tot_freq)

################################################################################
if __name__=='__main__':
    main()
