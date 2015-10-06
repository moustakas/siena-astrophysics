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
    comdist=cosmo.comoving_distance(redshift).value

    # Convert spherical to Cartesian Coordinates
    x=comdist*np.sin(dec)*np.cos(ra)
    y=comdist*np.sin(dec)*np.sin(ra)
    z=comdist*np.cos(dec)

    coords=np.column_stack((x,y,z))

    return coords


# This is the way we think we should calculate para and perp.

def our_para_perp(r0,r1):
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

    # First compute R_LOS and dR
    R_LOS = (r0 + r1)/2
    dR = r1 - r0
    R_LOS_mag = mag(R_LOS)

    # Dot product
    
    # We do this here ([:,0], e.g.) because we expect r1 to be an array.
    R_para = (dR[:,0]*R_LOS[:,0] + dR[:,1]*R_LOS[:,1] + dR[:,2]*R_LOS[:,2])/R_LOS_mag
    
    dR_mag = mag(dR)
    
    # Make use of the Pythagorean theorem
    R_perp = np.sqrt(dR_mag*dR_mag - R_para*R_para)
    
    #print i,lo1,indexlo,indexhi,len(R_para),len(paras)

    return R_para,R_perp


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



# This is for 1D.

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

    # Because we know that r1 is an array.
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
        ra=(data['PLUG_RA'])*((math.pi)/180)
        dec=((math.pi)/2)-((data['PLUG_DEC'])*((math.pi)/180))
        redshift = data['Z']

        #del data
    else:
        # Opening txt file (Mocks) 'b'
        print 'Reading in Text File'
        r=np.loadtxt(infilename)

        # Radians
        ra=(r[:,0])*((math.pi)/180)
        dec=((math.pi)/2)-((r[:,1])*((math.pi)/180))
        redshift=r[:,2]

        #del r

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

    if xyz:
        r=np.loadtxt(infilename)
        coords = np.column_stack((r[:,0],r[:,1],r[:,2]))
        del r,h1
    if return_radecz:
        coords = np.column_stack((np.rad2deg(ra),np.rad2deg(-(dec-np.pi/2)),redshift))

    else:
        coords = radecredshift2xyz(ra,dec,redshift)
        del ra,dec,redshift

    return coords





