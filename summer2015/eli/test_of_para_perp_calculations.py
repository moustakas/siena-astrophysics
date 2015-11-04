############ DISTANCE TEST ##################
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.io
import matplotlib.pylab as plt
import math
import matplotlib.pylab as plt

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Distance



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

    print "comdist: "
    print cosmo.comoving_distance(redshift)

    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, distance=comdist*u.Mpc)

    halfpi = math.pi/180.
    # Convert degrees to radians and spherical coordinates
    RArad=np.deg2rad(ra)
    #Decrad=((dec)*((math.pi)/180.))#((math.pi)/2.)-
    Decrad=halfpi-np.deg2rad(dec)
    #Decrad=np.deg2rad(dec)

    # Convert spherical to Cartesian Coordinates
    x=comdist*np.sin(Decrad)*np.cos(RArad)
    y=comdist*np.sin(Decrad)*np.sin(RArad)
    z=comdist*np.cos(Decrad)

    print np.sqrt(x*x + y*y + z*z)
    xc = c.icrs.cartesian.x
    yc = c.icrs.cartesian.y
    zc = c.icrs.cartesian.z
    print np.sqrt(xc*xc + yc*yc + zc*zc)

    print x[0],xc[0]
    print y[0],yc[0]
    print z[0],zc[0]

    coords=np.column_stack((x,y,z))
    #coords=np.column_stack((xc.value,yc.value,zc.value))

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
paras1 = []
perps1 = []
for i in range(0,ngals):
    for j in range(i+1,ngals):

        r1=coords[i]     
        r2=coords[j]

        r1mag = mag(r1)
        r2mag = mag(r2)
        rat = r1mag/r1mag

        #print "mags:",r1mag,r2mag

        # First compute R_LOS (line-of-sight vector) and dR
        R_LOS = (r1 + r2)/2.
        #print "DIFF: ",R_LOS-(r1 + r2*rat)
        #R_LOS = r1 + r2*rat
        #dR = r2 - r1
        dR = r1 - r2
        R_LOS_mag = mag(R_LOS)

        #print "dR:", dR
        #print "R_LOS_mag: ",R_LOS_mag

        # Dot product (normalized)        
        R_para = (dR[0]*R_LOS[0] + dR[1]*R_LOS[1] + dR[2]*R_LOS[2])/R_LOS_mag
        dR_mag = mag(dR) 
        #R_para=np.absolute(R_para)
        # Make use of the Pythagorean theorem        
        R_perp = np.sqrt(dR_mag*dR_mag - R_para*R_para)
                
        print "----------------"
        print "gal1: %f %f %f" % (r1[0],r1[1],r1[2])
        print "gal2: %f %f %f" % (r2[0],r2[1],r2[2])
        print "para/perp: %f %f" % (R_para,R_perp)
        paras1.append(R_para)
        perps1.append(R_perp)

###############################################################################
# Lado's Code
###############################################################################
paras2=[]
perps2=[]
for i in range(0,ngals):
    for j in range(i+1,ngals):
    
        #print "----------"
        x1=coords[i,0]
        x2=coords[j,0]
        y1=coords[i,1]
        y2=coords[j,1]
        z1=coords[i,2]
        z2=coords[j,2]

        #print x1,y1,z1
        #exit()

        d1 = x1-x2
        d2 = y1-y2
        d3 = z1-z2
        r2l = d1*d1 + d2*d2 + d3*d3

        #print "dR", d1,d2,d3

        gd1 = np.sqrt((x1**2)+(y1**2)+(z1**2))
        gd2 = np.sqrt((x2**2)+(y2**2)+(z2**2))
        rat = gd1/gd2

        #print "mags:",gd1,gd2

        xb = x1 + x2*rat
        yb = y1 + y2*rat
        zb = z1 + z2*rat

        db2 = xb*xb + yb*yb + zb*zb
        #print "dbmag: ",np.sqrt(db2)

        mu = np.absolute((xb*d1 + yb*d2 + zb*d3)/np.sqrt(r2l)/np.sqrt(db2))
        rr = np.sqrt(r2l)

        rpar=rr*mu
        rperp=rr*np.sqrt(1-(mu*mu))

        #print "rpar,rperp:",rpar,rperp
        
        paras2.append(rpar)
        perps2.append(rperp)    



paras1 = np.abs(np.array(paras1))
paras2 = np.array(paras2)

perps1 = np.array(perps1)
perps2 = np.array(perps2)


plt.figure(figsize=(16,6))
plt.subplot(1,2,1)
plt.plot(paras1,(paras1-paras2)/paras1,'o')
plt.xlabel('Our Code')
plt.ylabel("Fractional difference")
plt.title('Paras')


plt.subplot(1,2,2)
plt.plot(perps1,(perps1-perps2)/perps1,'o')
plt.xlabel('Our Code')
plt.ylabel("Fractional difference")
plt.title('Perps')

plt.tight_layout()

plt.show()



        
