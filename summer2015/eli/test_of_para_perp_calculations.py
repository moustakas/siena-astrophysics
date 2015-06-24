############ DISTANCE TEST ##################
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.io
import matplotlib.pylab as plt
import math

################################################################################
# Get magnitude of a vector
################################################################################
def mag(vec):

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

    # Convert degrees to radians and spherical coordinates
    RArad=(ra)*((math.pi)/180.)
    Decrad=((math.pi)/2.)-((dec)*((math.pi)/180.))

    # Convert spherical to Cartesian Coordinates
    x=comdist*np.sin(Decrad)*np.cos(RArad)
    y=comdist*np.sin(Decrad)*np.sin(RArad)
    z=comdist*np.cos(Decrad)

    coords=np.column_stack((x,y,z))

    return coords


################################################################################
# 10 Galaxies (ra,dec,z)    
################################################################################
galdat=np.array([[221.42511802,6.48006286,0.5891996],
                [143.26441531,16.56591479,0.68778771],
                [250.39513305,32.56592786,0.51523036],
                [192.75572856,-3.06978128,0.56942481],
                [136.94757689,0.61384908,0.53580755],
                [126.64582688,13.7770392,0.54264402],
                [127.16685388,42.06680514,0.53326231],
                [130.6626883,4.63926095,0.62787837],
                [242.14836333,12.94153412,0.46952653],
                [236.24144054,23.61780059,0.53417909]])

ra=galdat[:,0] # Right Ascension Values
dec=galdat[:,1] # Declination Values
z=galdat[:,2] # Redshift

coords = radecredshift2xyz(ra,dec,z)


################################################################################
# Do all the calculations!
#
# Parallel and perpendicular distances
################################################################################

ngals = len(coords)

for i in range(0,ngals):
    for j in range(i+1,ngals):

        r1=coords[i]     
        r2=coords[j]

        # First compute R_LOS (line-of-sight vector) and dR
        R_LOS = (r1 + r2)/2.
        dR = r2 - r1
        R_LOS_mag = mag(R_LOS)

        # Dot product (normalized)        
        R_para = (dR[0]*R_LOS[0] + dR[1]*R_LOS[1] + dR[2]*R_LOS[2])/R_LOS_mag
        dR_mag = mag(dR)

        # Make use of the Pythagorean theorem        
        R_perp = np.sqrt(dR_mag*dR_mag - R_para*R_para)
                
        print "----------------"
        print "gal1: %f %f %f" % (r1[0],r1[1],r1[2])
        print "gal2: %f %f %f" % (r2[0],r2[1],r2[2])
        print "para/perp: %f %f" % (R_para,R_perp)



