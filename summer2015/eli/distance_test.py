############ DISTANCE TEST ##################
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.io
import matplotlib.pylab as plt
import math

# 10 Galaxies (ra,dec,z)    
galdat=np.array([[221.42511802,6.48006286,0.5891996],[143.26441531,16.56591479,0.68778771],[250.39513305,32.56592786,0.51523036],[192.75572856,-3.06978128,0.56942481],[136.94757689,0.61384908,0.53580755],[126.64582688,13.7770392,0.54264402],[127.16685388,42.06680514,0.53326231],[130.6626883,4.63926095,0.62787837],[242.14836333,12.94153412,0.46952653],[236.24144054,23.61780059,0.53417909]])
ra=galdat[:,0] # Right Ascension Values
dec=galdat[:,1] # Declination Values
z=galdat[:,2] # Redshift

# Comoving Distances In Mpc
cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
comdist=cosmo.comoving_distance(z)

# Converting RA and Dec to Radians
RArad=(ra)*((math.pi)/180)
Decrad=((math.pi)/2)-((dec)*((math.pi)/180))

# Converting to Cartesian Coordinates
x=comdist*np.sin(Decrad)*np.cos(RArad)
y=comdist*np.sin(Decrad)*np.sin(RArad)
z=comdist*np.cos(Decrad)
coords=np.column_stack((x,y,z))

def mag(vec):

        m = None
        # First check if it is an 3xn array of coordinates....
        if type(vec[0])==np.ndarray or type(vec[0])==astropy.units.quantity.Quantity:
            m = np.sqrt(vec[:,0]**2 + vec[:,1]**2 + vec[:,2]**2)
        else:
            # Or if it is just the 3 coordinates.
            m = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)

        return m
    ################################################################################
ngals = len(coords)


paras = []
perps = []

for i in range(0,ngals):
    if i!=ngals-1:
        # First compute R_LOS and dR
        r1=coords[i]     
        R_LOS1 = (r1 + coords[i+1:])/2.
        dR1 = coords[i+1:] - r1
        R_LOS_mag1 = mag(R_LOS1)

        # Dot product        
        R_para1 = (dR1[:,0]*R_LOS1[:,0] + dR1[:,1]*R_LOS1[:,1] + dR1[:,2]*R_LOS1[:,2])/R_LOS_mag1  
        dR_mag1 = mag(dR1)

        # Make use of the Pythagorean theorem        
        R_perp1 = np.sqrt(dR_mag1*dR_mag1 - R_para1*R_para1)
               
                
        paras += R_para1.tolist()
                
        perps += R_perp1.tolist()



#### Final distances are contained in the lists named paras and perps. #### 



