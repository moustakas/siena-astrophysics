# location.py -- some nice location handling utilities
# [System]

# [Installed]
import numpy as np
import pywcs
import cosmolopy
import scipy.interpolate

# [Package]

# [Constants]

def _cosmology():
    '''returns the values that match kCorrect and idl'''
    p = cosmolopy.fidcosmo # default values from cosmolpy
    n  = {
        'h': 1.0,
        'omega_M_0':0.27, 
        'omega_lambda_0':0.73,
        'sigma_8':0.8,
    }
    p.update(n)
    return p

def radec2xy(header, ra, dec):
    '''From a header, and radec position, return the x,y position relative to an image'''
    wcs = pywcs.WCS(header)
    x, y = wcs.wcs_sky2pix(ra,dec,0)
    return x,y
    # try:
    #     import pywcsgrid2 # slow so leave it here.
    #     t = pywcsgrid2.wcs_transforms.WcsSky2PixelTransform(header)
    #     ll = np.vstack((ra,dec)).T    
    #     xy = t.transform(ll)
    #     return xy[:,0], xy[:,1]
    # except 

def xy2radec(header, xx, yy=None):
    '''Returns (ra, dec) for a set of x and y values that match a wcs header.
    you can pass in xy = np.zeros((n,2)) array as xx and leave yy as None'''
    wcs = pywcs.WCS(header)
    if yy is None:
        radec = wcs.wcs_pix2sky(xx,0)
        return radec[:,0], radec[:,1]
    else:
        return wcs.wcs_pix2sky(xx,yy,0)
    
    # if yy is None:
    #     xy = xx
    # else:
    #     xy = np.vstack((xx,yy)).T
    # import pywcsgrid2 # slow so leave it here
    # t = pywcsgrid2.wcs_transforms.WcsPixel2SkyTransform(header)
    # radec = t.transform(xy)
    # return radec[:,0], radec[:,1]



def convert2distance(ra, dec, z, center):
    '''Convert ra,dec,z into Mpc/h units
    
    '''
    dra = np.radians(ra-center[0])
    ddec = np.radians(dec-center[1])
    
    p = _cosmology()
    dt = cosmolopy.distance.comoving_distance_transverse(z, **p)*p['h']
    rx = dt*dra*np.cos(np.radians(ddec))
    ry = dt*ddec
    rz = cosmolopy.distance.comoving_distance(z, **p)*p['h']
    return rx,ry,rz
    # p = _cosmology()
    # dt = cosmolopy.distance.comoving_distance_transverse(z, **p)*p['h']
    # rx = dt*np.radians(ra - center[0])*np.cos(np.radians(center[1]))
    # ry = dt*np.radians(dec - center[1])
    # rz = cosmolopy.distance.comoving_distance(z, **p)*p['h']
    # return rx,ry,rz

def comoving_distance(z):
    p = _cosmology()
    rz = cosmolopy.distance.comoving_distance(z, **p)*p['h']
    return rz

def redshift_interpol():
    '''return a function that takes rz comoving distances and converts
    them back to nice redshifts'''
    z = np.arange(0,5,0.01)
    rz = comoving_distance(z)
    f = scipy.interpolate.interp1d(rz,z)
    return f


def redshiftvolume(area, zrange, zbins=1000):
    '''Gets the volume by integrating the area times zbins
    area is in deg^2 
    zrange = [zmin, zmax] 
    zbins = number of redshift bins
    returns the volume in (Mpc/h)^3
                 /-------------+
        /-------               |
    +--                        |
    |                          | area
    +--                        |
        \-------               |
                 \-------------+
    zmin                    zmax
    '''
    
    # setup the variables
    p = _cosmology()
    z = np.linspace(zrange[0], zrange[1], zbins)
    z2 = np.linspace(zrange[0], zrange[1], zbins+1)
    
    
    dt = cosmolopy.distance.comoving_distance_transverse(z, **p)*p['h']
    
    area = (dt*np.radians(np.sqrt(area)))**2
    dz = np.diff(cosmolopy.distance.comoving_distance(z2, **p)*p['h'])
    volume = np.sum(area*dz)
    return volume


