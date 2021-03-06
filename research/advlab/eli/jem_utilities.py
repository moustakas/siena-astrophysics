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

# Test of repo

def mag(vec):
    """Get magnitude of a vector.

    Args:
      vec (numpy.ndarray): X,Y, and Z components of a vector

    Returns:
      m (numpy.ndarray): The magnitude of the vector.

    """

    m = None
    if type(vec)==np.ndarray:
        if vec.shape==(3,):
            m = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
        else:
            m = np.sqrt(vec[:,0]**2 + vec[:,1]**2 + vec[:,2]**2)
    else:
        m = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

    return m




# Converting RA and Dec and redshift to Cartesian

def radecredshift2xyz(ra,dec,redshift):
    """Convert arrays of Right Ascension, Declination, and
                        Redshift to Cartesian Coordinates.

    Args:
        ra (numpy.ndarray): Right Ascension Values
        dec (numpy.ndarray): Declination Values
        redshift (numpy.ndarray): Redshift Values

    Returns:
        coords (numpy.ndarray): A three column array with
                                correspoding values of X,Y,
                                and Z
    """

    # Comoving Distances In Mpc
    cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
    comdist=cosmo.comoving_distance(redshift).value * 0.7 # Trying 0.7 for Lado's code.

    # Convert spherical to Cartesian Coordinates
    #x=comdist*np.sin(dec)*np.cos(ra)
    #y=comdist*np.sin(dec)*np.sin(ra)
    #z=comdist*np.cos(dec)

    # Reproducing Lado's stuff.
    x=comdist*np.cos(dec)*np.cos(ra)
    y=comdist*np.cos(dec)*np.sin(ra)
    z=comdist*np.sin(dec)

    coords=np.column_stack((x,y,z))

    return coords

################################################################################
################################################################################
def radecredshift2xyz_with_weights(oldcoords):

    ra = oldcoords[:,0]
    dec = oldcoords[:,1]
    redshift = oldcoords[:,2]
    weights = oldcoords[:,3]

    # Comoving Distances In Mpc
    cosmo=FlatLambdaCDM(H0=70,Om0=0.274)
    comdist=cosmo.comoving_distance(redshift).value * 0.7 # Trying 0.7 for Lado's code.

    # Convert spherical to Cartesian Coordinates
    #x=comdist*np.sin(dec)*np.cos(ra)
    #y=comdist*np.sin(dec)*np.sin(ra)
    #z=comdist*np.cos(dec)

    # Reproducing Lado's stuff.
    x=comdist*np.cos(dec)*np.cos(ra)
    y=comdist*np.cos(dec)*np.sin(ra)
    z=comdist*np.sin(dec)

    coords=np.column_stack((x,y,z,weights,comdist))

    return coords

################################################################################
################################################################################


# This is the way we think we should calculate para and perp.

def our_para_perp_with_weights(r1,r2):
    """Calculate r_parallel and r_perpendicular distances
                                    between two galaxies.

    Args:
        r0 (numpy.ndarray): A three dimensional vector to
                            galaxy 1
        r1 (numpy.ndarray): A three dimensional vector to
                            galaxy 2

    Returns:
        R_para () : The magnitude of the r_parallel distance
        R_perp () : The magnitude of the r_perpendicular distance
    """
    w1 = r1[3]
    w2 = r2[:,3]
    # First compute R_LOS and dR
    R_LOS = (r1 + r2)/2
    dR = r2 - r1
    R_LOS_mag = mag(R_LOS)

    # Dot product
    
    # We do this here ([:,0], e.g.) because we expect r1 to be an array.
    R_para = (dR[:,0]*R_LOS[:,0] + dR[:,1]*R_LOS[:,1] + dR[:,2]*R_LOS[:,2])/R_LOS_mag
    
    dR_mag = mag(dR)
    
    # Make use of the Pythagorean theorem
    R_perp = np.sqrt(dR_mag*dR_mag - R_para*R_para)
    
    #print i,lo1,indexlo,indexhi,len(R_para),len(paras)
    weights = w1*w2
    return R_para,R_perp,weights


# This is the way we think Lado calculates para and perp.

def lado_para_perp(r1,r2):
    """Calculate r_parallel and r_perpendicular distances
                        between two galaxies, Lado's way.

    Args:
        r1 (numpy.ndarray): A three dimensional vector to
                                                 galaxy 1
        r2 (numpy.ndarray): A three dimensional vector to
                                                 galaxy 2

    Returns:
        rpar () : The magnitude of the r_parallel distance
        rperp () : The magnitude of the r_perpendicular distance
    """

    #x1=r1[:,0]
    #y1=r1[:,1]
    #z1=r1[:,2]
    x1=r1[0]
    y1=r1[1]
    z1=r1[2]

    # Because we know that r1 is an array.
    x2=r2[:,0]
    y2=r2[:,1]
    z2=r2[:,2]

    d1 = x1-x2
    d2 = y1-y2 
    d3 = z1-z2
    
    r2l = d1*d1 + d2*d2 + d3*d3

    gd1 = mag(r1)
    gd2 = mag(r2)
    rat = gd1/gd2

    xb = x1 + (x2)*rat
    yb = y1 + (y2)*rat
    zb = z1 + (z2)*rat

    db2 = xb*xb + yb*yb + zb*zb

    mu = np.absolute(((xb*d1 + yb*d2 + zb*d3)/np.sqrt(r2l)/np.sqrt(db2)))
    rr = np.sqrt(r2l)
     
    rpar=rr*mu
    rperp=rr*np.sqrt(1-(mu*mu))

    return rpar,rperp



################################################################################
def angular_sep(r1,d1,r2,d2):
    # Assumes ra and dec are in radians
    #cos(A) = sin(d1)sin(d2) + cos(d1)cos(d2)cos(ra1-ra2)

    cosA = np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(r1-r2)

    A = np.arccos(cosA)

    return A

################################################################################
# This is for 1D new way of calcing distance. Discussion with Rose?
################################################################################
def one_dimension_trial(r1,r2):
    """Calculates standard distance between two galaxies 

    Args:
        r1 (numpy.ndarray): A three dimensional vector to
                                                 galaxy 1
        r2 (numpy.ndarray): A three dimensional vector to
                                                 galaxy 2

    Returns:
        distances (): The standard distance between two galaxies
        fake_vals (numpy.ndarray): An array of zeros, this
                                   is needed to build an array
                                   shape.  
    """
    cosmo=FlatLambdaCDM(H0=70,Om0=0.3)

    ra1=r1[0]
    dec1=r1[1]
    z1=r1[2]
    d1=r1[3]

    # Because we know that r1 is an array.
    ra2=r2[:,0]
    dec2=r2[:,1]
    z2=r2[:,2]
    d2=r2[:,3]

    asep = angular_sep(ra1,dec1,ra2,dec2)
    asep_in_min = np.rad2deg(asep)*60.

    #avgz = (z1+z2)/2.

    x = cosmo.kpc_comoving_per_arcmin(z1).value
    #x = cosmo.kpc_comoving_per_arcmin(avgz).value
    x *= asep_in_min/1000.0 # Convert to Mpc

    #d1 = cosmo.comoving_distance(z1).value
    #d2 = cosmo.comoving_distance(z2).value

    y = d2-d1

    distances = np.sqrt(x*x + y*y)

    fake_vals = np.zeros(len(distances))

    return distances,fake_vals

################################################################################

    
# This is for 1D.

################################################################################
################################################################################
def one_dimension_with_weights(r1,r2):
    x1=r1[0]
    y1=r1[1]
    z1=r1[2]
    w1=r1[3]

    # Because we know that r2 is an array.
    x2=r2[:,0]
    y2=r2[:,1]
    z2=r2[:,2]
    w2=r2[:,3]

    d1 = x1-x2
    d2 = y1-y2 
    d3 = z1-z2

    distances = mag([d1,d2,d3])

    weights = w1*w2

    return distances,weights
################################################################################
    

def one_dimension(r1,r2):
    """Calculates standard distance between two galaxies 

    Args:
        r1 (numpy.ndarray): A three dimensional vector to
                                                 galaxy 1
        r2 (numpy.ndarray): A three dimensional vector to
                                                 galaxy 2

    Returns:
        distances (): The standard distance between two galaxies
        fake_vals (numpy.ndarray): An array of zeros, this
                                   is needed to build an array
                                   shape.  
    """
    x1=r1[0]
    y1=r1[1]
    z1=r1[2]

    # Because we know that r2 is an array.
    x2=r2[:,0]
    y2=r2[:,1]
    z2=r2[:,2]

    d1 = x1-x2
    d2 = y1-y2 
    d3 = z1-z2

    distances = mag([d1,d2,d3])

    fake_vals = np.zeros(len(distances))

    return distances,fake_vals
    

# Using PySurvey

def pysurvey_distance(r1,r2):
    """Calculates standard distance between two
                        galaxies using pysurvey

    Args:
        r1 (numpy.ndarray): A three dimensional
                            vector to galaxy 1
        r2 (numpy.ndarray): A three dimensional
                            vector to galaxy 2

    Returns:
        dist (): The standard distance between
                                  two galaxies
        fake_vals (numpy.ndarray): An array of zeros,
                                   this is needed to
                                   build an array shape.  
    """
    ra1=r1[0]
    dec1=r1[1]
    z1=r1[2]

    ra2=r2[:,0]
    dec2=r2[:,1]
    z2=r2[:,2]

    loc1 = location.convert2distance(ra1,dec1,z1,[0.,0.])
    loc2 = location.convert2distance(ra2,dec2,z2,[0.,0.])

    dist = mag([loc1[0]-loc2[0],loc1[1]-loc2[1],loc1[2]-loc2[2]])

    fake_vals = np.zeros(len(dist))

    return  dist,fake_vals


def get_coordinates(infilename,xyz=False,maxgals=0,return_radecz=False):
    """Grabs either X.Y,Z coordinates or RA,DEC,Z values
                                          of a data set.

    Args:
        infilename: The name of the data file.
        xyz (Boolean): If true, the input file is one
                         with xyz coordinates already 
        maxgals (int): The number of galaxies to be
                       calculated, if 0, then all 
                       galaxies will be calculated.
        return_radecz (): If true, Right Ascension,
                          Declination, and Redshift
                          will be returned.
    Returns:
        coords (numpy.ndarray): Either an array of
                                   XYZ or RA,DEC,Z
    """
    isdatafile = False
    if(infilename.find('fits')>=0):
        isdatafile = True

    if isdatafile:
        # Opening FITS file (SDSS Data) 'a' 
        print 'Reading in FITS Data'
        hdulist1=fits.open(infilename)
        hdulist1.info()
        h1=hdulist1[1]
        data=h1.data
        del h1


        # Radians
        #ra=(data['PLUG_RA'])*((math.pi)/180)
        #dec=((math.pi)/2)-((data['PLUG_DEC'])*((math.pi)/180))
        ra=np.deg2rad(data['PLUG_RA'])
        dec=np.deg2rad(data['PLUG_DEC'])
        redshift = data['Z']

        del data
    else:
        # Opening txt file (Mocks) 'b'
        print 'Reading in Text File'
        r=np.loadtxt(infilename)
        if xyz:
            coords = np.column_stack((r[:,0],r[:,1],r[:,2]))
            del r
            return coords
            exit()
        else:    
            # Radians
            #ra=(r[:,0])*((math.pi)/180)
            #dec=((math.pi)/2)-((r[:,1])*((math.pi)/180))
            ra=np.deg2rad(r[:,0])
            dec=np.deg2rad(r[:,1])
            redshift=r[:,2]

            del r

    # Made some common cuts
    index0 = redshift<0.7
    index1 = redshift>0.43
    index = index0*index1

    '''
    # Fiducial cuts?
    indexcut = dec<0.35
    indexcut *= ra<3.7
    indexcut += ra>3.7
    index *= indexcut
    '''

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

    if xyz:
        r=np.loadtxt(infilename)
        coords = np.column_stack((r[:,0],r[:,1],r[:,2]))
        del r
    if return_radecz:
        #coords = np.column_stack((np.rad2deg(ra),np.rad2deg(-(dec-np.pi/2)),redshift))
        cosmo=FlatLambdaCDM(H0=70,Om0=0.310) # Doing 0.310 to match Lado's code.
        d = cosmo.comoving_distance(redshift).value * 0.7 # Trying the 0.7 to match Lado's code.
        #coords = np.column_stack((np.deg2rad(ra),np.deg2rad(dec),redshift,d))
        #print ra[0:100]
        #print dec[0:100]
        coords = np.column_stack((ra,dec,redshift,d))

    else:
        coords = radecredshift2xyz(ra,dec,redshift)
        del ra,dec,redshift

    return coords

################################################################################
################################################################################
def get_coordinates_with_weight(infilename):

    r=np.loadtxt(infilename,unpack=True)

    ra=np.deg2rad(r[0])
    dec=np.deg2rad(r[1])
    redshift=r[2]
    #weights=r[3]
    #weights=r[4]
    weights=np.ones(len(redshift))
    #weights=0.1*np.ones(len(redshift))

    del r

    # Made some common cuts
    index0 = redshift<0.7
    index1 = redshift>0.43
    index = index0*index1

    ra = ra[index]
    dec = dec[index]
    redshift = redshift[index]
    weights = weights[index]
    coords = np.column_stack((ra,dec,redshift,weights))

    return coords

################################################################################
################################################################################
                                                                                
def corr_est(DD,DR,RR,ngals,nrands,nbins,
             oned=False):
    """Calculates the Landy-Szalay Correlation Function Estimator,
                                     given three frequency arrays

        Args:
            DD (numpy.ndarray or ASCII file) : The frequency array
                                               of data-data distances
            DR (numpy.ndarray or ASCII file) : The frequency array of
                                                data-random distances
            RR (numpy.ndarray or ASCII file) : The frequency array of
                                               random-random distances
            ngals (int) : The number of data points originally
                                                     utilized
            nrands (int) : The number of random points originally
                                                         utilized
        Returns:
            Xi (numpy.ndarray) : The frequency array of the Correlation
                                                              Estimator
    """
    
    if type(DD)==str:
        D_D=np.loadtxt(DD,dtype='float')
    else:
        D_D=DD

    if type(DR)==str:
        D_R=np.loadtxt(DR,dtype='float')
    else:
        D_R=DR

    if type(RR)==str:
        R_R=np.loadtxt(RR,dtype='float')
    else:
        R_R=RR

    D_D = D_D.transpose()
    R_R = R_R.transpose()
    D_R = D_R.transpose()

    D_D+=np.fliplr(D_D)
    D_R+=np.fliplr(D_R)
    R_R+=np.fliplr(R_R)
    
    D_D /=(ngals**2-ngals)/2.
    D_R /=(nrands*ngals)/1.
    R_R /=(nrands**2-nrands)/2.
    Xi = ((D_D - 2*D_R + R_R)/(R_R+(R_R==0)))*(R_R!=0)

    if oned:
        index=nbins/2
        Xi1 = Xi.transpose()[index]
        Xi = Xi1[index:]
    print Xi, len(Xi)        
    return Xi

def corr_plot(Xi,x0,x1,y0,y1,title,xlab,ylab,oned=False):
    """Provides histogram and scatter plots for a number of different
                                                             infiles.
        Args:
            infile (numpy.ndarry) : The data to be plotted
            x0 (int) : Lower x bound
            x1 (int) : Upper x bound
            y0 (int) : Lower y bound
            y1 (int) : Upper y bound
            title (str) : The title of the plot
            xlab (str) : Label for the x-axis
            ylab (str) : Label for the y axis
            1D (Boolean) : Set to true if 1D Correlation plot is desired
        Returns:
            A plot of the data

    """
    import numpy as np
    import matplotlib.pylab as plt
        
    if oned==True:
        npts=len(Xi)
        xvals = np.linspace(0,x1,npts)
        fig= plt.figure()
        plt.plot(xvals,Xi,'o')
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        plt.savefig('onedcorr.png')
        
    else:
        Xi+=np.fliplr(Xi)
        fig=plt.figure(figsize=(8,8))
        extent=[x0,x1,y0,y1]
        plot=plt.imshow(Xi,extent=extent)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        plt.show()
        

def twopoint_hist(infile1,infile2,nbins,rangeval,
                         range1=None,range2=None,
                                     oned=False):
    """Given command-line arguments, will return
        frequency arrays for galactic distances.
    Args:
        See Below
    Returns:
        Dependant on command-line arguments.
    """
    infilename0 = infile1
    infilename1 = infile2

    range1lo = None
    range1hi = None
    if range1 is not None:
        range1lo = int(range1.split('-')[0])
        range1hi = int(range1.split('-')[1])

    range2lo = None
    range2hi = None
    if range2 is not None:
        range2lo = int(range2.split('-')[0])
        range2hi = int(range2.split('-')[1])

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
    nbins=nbins
    rangeval=rangeval

    tot_freq = np.zeros((nbins,nbins)) 

    # Figure out the chunking.
    

    chunk_size = 50
    nchunks = len(coords0cut)/chunk_size     #ngals_for_calculation/chunk_size

    ncalcs_per_chunk = chunk_size*len(coords1cut) #chunk_size*ngals1

    # These will store the calculations.
    paras = np.zeros(ncalcs_per_chunk)
    perps = np.zeros(ncalcs_per_chunk)

    indexlo = 0
    indexhi = 0
    
    if samefile and (infilename0.find('fits')>=0):
        calc = "DD"
    elif samefile and (infilename0.find('fits')==-1):
        calc = "RR"
    else:
        calc = "DR"
        
    #Calculation Loop
    for j in xrange(nchunks):
        lo = j*chunk_size
        hi = (j+1)*chunk_size
        print "Performing calculations for %s %d chunk: %d-%d" % (calc,j,lo,hi)

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

            if oned:
                temp_paras,temp_perps = one_dimension(r0,other_gals)
            else:
                
                temp_paras,temp_perps = our_para_perp(r0,other_gals)

            paras[indexlo:indexhi] = temp_paras
            perps[indexlo:indexhi] = temp_perps

            indexlo = indexhi
        
        # Histogram the values.
        hist=plt.hist2d(perps[0:indexhi],paras[0:indexhi],bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
        tot_freq += hist[0]

        tot_freq += hist[0]

        indexlo=0
        indexhi=0

        del hist

    print tot_freq
    return tot_freq

    

################################################################################
def define_ranges(loranges, hiranges, maxsep=200):

    ncoords = len(loranges)
    ngrids = []
    gridwidths = []

    #print "loranges and hiranges"
    #print loranges
    #print hiranges

    for i in range(0,ncoords):

        r = hiranges[i]-loranges[i];

        ngrids.append(int(r/maxsep))
        gridwidths.append(r/ngrids[i])

    return ngrids,gridwidths


################################################################################
################################################################################
def assign_grid_coordinate(coords, loranges, hiranges, gridwidths):

    # Number of points is number of x-axes.
    ncoords = len(coords[:,0])

    grid_coordinates = -1*np.ones((3,ncoords),dtype=int)

    for i in range(0,3):
        grid_coordinates[i] = ((coords[:,i]-loranges[i])/gridwidths[i]).astype(int)
        # Look for the ones on the edge and move them down one.
        index = (coords[:,i]==hiranges[i])
        grid_coordinates[i][index] -= 1
        #print grid_coordinates[i][index]
        #exit()

    # Need to check for anything less than 0!!!!

    return grid_coordinates

################################################################################
################################################################################

def twopoint_hist_grid(infile1,infile2,nbins,rangeval,
                         range1=None,range2=None,
                                     oned=False):
    """Given command-line arguments, will return
        frequency arrays for galactic distances.

    Args:
        See Below

    Returns:
        Dependant on command-line arguments.
    """

    infilename0 = infile1
    infilename1 = infile2

    range1lo = None
    range1hi = None
    if range1 is not None:
        range1lo = int(range1.split('-')[0])
        range1hi = int(range1.split('-')[1])

    range2lo = None
    range2hi = None
    if range2 is not None:
        range2lo = int(range2.split('-')[0])
        range2hi = int(range2.split('-')[1])

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

    coords0 = get_coordinates(infilename0,False,0,return_radecz=True)
    coords1 = get_coordinates(infilename1,False,0,return_radecz=True)

    # Convert to x,y,z
    print "Converting to x,y,z...."
    coords0 = radecredshift2xyz(coords0[:,0],coords0[:,1],coords0[:,2])
    coords1 = radecredshift2xyz(coords1[:,0],coords1[:,1],coords1[:,2])

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
   
    ############################################################################
    # Break up into voxels.
    ############################################################################
    loranges = [np.min((coords0cut[:,0],coords1cut[:,0])),np.min((coords0cut[:,1],coords1cut[:,1])),np.min((coords0cut[:,2],coords1cut[:,2]))]
    hiranges = [np.max((coords0cut[:,0],coords1cut[:,0])),np.max((coords0cut[:,1],coords1cut[:,1])),np.max((coords0cut[:,2],coords1cut[:,2]))]
    ngrids,gridwidths = define_ranges(loranges,hiranges, maxsep=200)

    grid_coords0 = assign_grid_coordinate(coords0cut, loranges, hiranges, gridwidths)
    grid_coords1 = assign_grid_coordinate(coords1cut, loranges, hiranges, gridwidths)

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
    #import pdb ; pdb.set_trace()
    # Fill them
    for i in range(0,ngals0):
        a,b,c = grid_coords0[0][i],grid_coords0[1][i],grid_coords0[2][i]
        voxels0[a][b][c].append(coords0[i])
    for i in range(0,ngals1):
        a,b,c = grid_coords1[0][i],grid_coords1[1][i],grid_coords1[2][i] 
        #print(a,b,c,i)
        #if i==537:
        #    import pdb ; pdb.set_trace()
        voxels1[a][b][c].append(coords1[i])
    # Convert lists to arrays
    for i in range(0,ngrids[0]):
        for j in range(0,ngrids[1]):
            for k in range(0,ngrids[2]):
                voxels0[i][j][k] = np.array(voxels0[i][j][k])

    for i in range(0,ngrids[0]):
        for j in range(0,ngrids[1]):
            for k in range(0,ngrids[2]):
                voxels1[i][j][k] = np.array(voxels1[i][j][k])


    ############################################################################
    # This is for the histogram.
    nbins=nbins
    rangeval=rangeval

    if oned==True:
        nbins*=2

    tot_freq = np.zeros((nbins,nbins)) 
    if oned==True:
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

                            if len(c0)>0 and len(c1)>0:
                                print ii,jj,kk, aa, bb, cc, len(c0),len(c1)

                            if len(c0)==0 or len(c1)==0:
                                continue

                            # These will store the calculations.
                            for index,r0 in enumerate(c0):

                                if samefile and aa==ii and bb==jj and cc==kk:
                                    paras,perps = one_dimension(r0,c1[index+1:])
                                else:
                                    paras,perps = one_dimension(r0,c1)

                                hist=np.histogram(paras,bins=nbins,range=(0,rangeval))

                                tot_freq += hist[0]

                                del hist

                            print tot_freq.sum()
   
    return tot_freq
    
################################################################################
################################################################################
def voxelize_the_data(coords0,coords1,maxsep=200):

    ngals0 = len(coords0[:,0])
    ngals1 = len(coords1[:,0])

    ############################################################################
    # Break up into voxels.
    ############################################################################
    loranges = [np.min((np.min(coords0[:,0]),np.min(coords1[:,0]))),np.min((np.min(coords0[:,1]),np.min(coords1[:,1]))),np.min((np.min(coords0[:,2]),np.min(coords1[:,2])))]
    hiranges = [np.max((np.max(coords0[:,0]),np.max(coords1[:,0]))),np.max((np.max(coords0[:,1]),np.max(coords1[:,1]))),np.max((np.max(coords0[:,2]),np.max(coords1[:,2])))]
    ngrids,gridwidths = define_ranges(loranges,hiranges, maxsep=maxsep)

    print "ngrids and gridwidths"
    print ngrids
    print gridwidths

    grid_coords0 = assign_grid_coordinate(coords0, loranges, hiranges, gridwidths)
    grid_coords1 = assign_grid_coordinate(coords1, loranges, hiranges, gridwidths)

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
        voxels0[a][b][c].append(coords0[i])
    for i in range(0,ngals1):
        a,b,c = grid_coords1[0][i],grid_coords1[1][i],grid_coords1[2][i] 
        voxels1[a][b][c].append(coords1[i])

    # Convert lists to arrays
    for i in range(0,ngrids[0]):
        for j in range(0,ngrids[1]):
            for k in range(0,ngrids[2]):
                voxels0[i][j][k] = np.array(voxels0[i][j][k])

    for i in range(0,ngrids[0]):
        for j in range(0,ngrids[1]):
            for k in range(0,ngrids[2]):
                voxels1[i][j][k] = np.array(voxels1[i][j][k])


    return voxels0,voxels1,ngrids,gridwidths,loranges,hiranges


############################################################################
def do_pair_counts(voxels0,voxels1,ngrids,nbins=10,maxrange=200,samefile=True):

    start_time_pc = time.time()

    voxel_combinations_so_far = []

    tot_distances = []

    tot_freq = np.zeros(nbins)

    tot_weight_val = 0

    tot_points_looped_over = 0

    ncalcs = 0

    weight_total = 0

    #Calculation Loop
    for ii in range(0,ngrids[0]):
        print ii,time.time()-start_time_pc
        for jj in range(0,ngrids[1]):
            for kk in range(0,ngrids[2]):

                c0 = voxels0[ii][jj][kk] 

                if len(c0)==0:
                    continue

                tot_points_looped_over += len(c0)

                #print 
                #print ii,jj,kk

                iimin = ii-1
                iimax = ii+2
                if iimax>=ngrids[0]:
                    iimax = ngrids[0]
                if iimin<0:
                    iimin=0
                jjmin = jj-1
                jjmax = jj+2
                if jjmax>=ngrids[1]:
                    jjmax = ngrids[1]
                if jjmin<0:
                    jjmin=0
                kkmin = kk-1
                kkmax = kk+2
                if kkmax>=ngrids[2]:
                    kkmax = ngrids[2]
                if kkmin<0:
                    kkmin=0

                for aa in range(iimin,iimax):
                    for bb in range(jjmin,jjmax):
                        for cc in range(kkmin,kkmax):

                            if samefile:
                                combination0 = "%03s%03s%03s" % (cc,bb,aa)
                                combination1 = "%03s%03s%03s" % (kk,jj,ii)
                                combination = "%s%s" % (max(combination0,combination1),min(combination0,combination1))

                                if combination in voxel_combinations_so_far:
                                    continue
                                else:
                                    voxel_combinations_so_far.append(combination)

                            c1 = voxels1[aa][bb][cc] 

                            if len(c1)==0:
                                continue

                            # These will store the calculations.
                            for index,r0 in enumerate(c0):

                                '''
                                if samefile and aa==ii and bb==jj and cc==kk:
                                    distances = distances=scipy.spatial.distance.pdist(c0)
                                else:
                                    distances=scipy.spatial.distance.cdist(c0,c1)

                                weights = np.ones(len(distances))
                                '''

                                if samefile and aa==ii and bb==jj and cc==kk:
                                    distances,weights = one_dimension_with_weights(r0,c1[index+1:])
                                else:
                                    distances,weights = one_dimension_with_weights(r0,c1)

                                #print weights
                                hist=np.histogram(distances,weights=weights,bins=nbins,range=(0,maxrange))

                                #weight_total += weights.sum()

                                tot_freq += hist[0]

                                #print tot_freq

                                ncalcs += len(distances)
                                
                                weight_tot_one = np.sum(weights)

                                #print weight_tot_one 
                                #print weights

                                #print '..........'
                                #print 'tot'
                                tot_weight_val += weight_tot_one
                                #print "weight one"
                                #print weight_tot_one, tot_weight_val
                                del hist

                    #print tot_weight_val
   
    #print np.sort(tot_distances)
    #print "tot points looped over: %d" % (tot_points_looped_over)
    print "Total weights: %f" % (tot_weight_val)
    print "# calcs      : %f" % (ncalcs)
    #new_tot_freq = tot_freq/tot_weight_val
    new_tot_freq = tot_freq
    #new_tot_freq *= ncalcs
    #print new_tot_freq
    return new_tot_freq
#####################################################################################
def do_pair_counts_2d(voxels0,voxels1,ngrids,nbins=10,maxrange=200,samefile=True):

    start_time_pc = time.time()

    voxel_combinations_so_far = []

    tot_distances = []

    tot_freq = np.zeros(nbins)

    tot_weight_val = 0

    tot_points_looped_over = 0

    ncalcs = 0

    weight_total = 0

    #Calculation Loop
    for ii in range(0,ngrids[0]):
        print ii,time.time()-start_time_pc
        for jj in range(0,ngrids[1]):
            for kk in range(0,ngrids[2]):

                c0 = voxels0[ii][jj][kk] 

                if len(c0)==0:
                    continue

                tot_points_looped_over += len(c0)

                #print 
                #print ii,jj,kk

                iimin = ii-1
                iimax = ii+2
                if iimax>=ngrids[0]:
                    iimax = ngrids[0]
                if iimin<0:
                    iimin=0
                jjmin = jj-1
                jjmax = jj+2
                if jjmax>=ngrids[1]:
                    jjmax = ngrids[1]
                if jjmin<0:
                    jjmin=0
                kkmin = kk-1
                kkmax = kk+2
                if kkmax>=ngrids[2]:
                    kkmax = ngrids[2]
                if kkmin<0:
                    kkmin=0

                for aa in range(iimin,iimax):
                    for bb in range(jjmin,jjmax):
                        for cc in range(kkmin,kkmax):

                            if samefile:
                                combination0 = "%03s%03s%03s" % (cc,bb,aa)
                                combination1 = "%03s%03s%03s" % (kk,jj,ii)
                                combination = "%s%s" % (max(combination0,combination1),min(combination0,combination1))

                                if combination in voxel_combinations_so_far:
                                    continue
                                else:
                                    voxel_combinations_so_far.append(combination)

                            c1 = voxels1[aa][bb][cc] 

                            if len(c1)==0:
                                continue

                            # These will store the calculations.
                            for index,r0 in enumerate(c0):

                                '''
                                if samefile and aa==ii and bb==jj and cc==kk:
                                    distances = distances=scipy.spatial.distance.pdist(c0)
                                else:
                                    distances=scipy.spatial.distance.cdist(c0,c1)

                                weights = np.ones(len(distances))
                                '''

                                if samefile and aa==ii and bb==jj and cc==kk:
                                    paras,perps,weights = our_para_perp_with_weights(r0,c1[index+1:])
                                else:
                                    paras,perps,weights = our_para_perp_with_weights(r0,c1)

                                #print weights
                                hist=np.histogram2d(perps,paras,bins=nbins,range=((-maxrange,maxrange),(-maxrange,maxrange)))

                                #weight_total += weights.sum()

                                tot_freq += hist[0]

                                #print tot_freq

                                ncalcs += len(distances)
                                
                                weight_tot_one = np.sum(weights)

                                #print weight_tot_one 
                                #print weights

                                #print '..........'
                                #print 'tot'
                                tot_weight_val += weight_tot_one
                                #print "weight one"
                                #print weight_tot_one, tot_weight_val
                                del hist

                    #print tot_weight_val
   
    #print np.sort(tot_distances)
    #print "tot points looped over: %d" % (tot_points_looped_over)
    print "Total weights: %f" % (tot_weight_val)
    print "# calcs      : %f" % (ncalcs)
    new_tot_freq = tot_freq/tot_weight_val
    #new_tot_freq = tot_freq
    #new_tot_freq *= ncalcs
    #print new_tot_freq
    return new_tot_freq

    

