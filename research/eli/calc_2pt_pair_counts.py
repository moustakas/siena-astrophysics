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


def get_coordinates(infilename,maxgals=0,return_radecz=False):
    """Grabs either X.Y,Z coordinates or RA,DEC,Z values
                                          of a data set.

    Args:
        infilename: The name of the data file.
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

        del r

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


    if return_radecz:
        coords = np.column_stack((np.rad2deg(ra),np.rad2deg(-(dec-np.pi/2)),redshift))

    else:
        coords = radecredshift2xyz(ra,dec,redshift)
        del ra,dec,redshift

    return coords






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

    coords0 = get_coordinates(infilename0,ngals_for_calculation,args.pysurvey)
    coords1 = get_coordinates(infilename1,nrands,args.pysurvey)
    
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
    nbins=200
    rangeval=200

    if args.oned:
        nbins*=2

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

    #Calculation Loop
    for j in xrange(nchunks):
        lo = j*chunk_size
        hi = (j+1)*chunk_size
        #print "Performing calculations for DR %d chunk: %d-%d" % (j,lo,hi)

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

            # Calc 1D using the pysurvey distance calc.
            if args.lado==False and args.oned==False and args.pysurvey==True:
                if len(other_gals)>0:
                    temp_paras,temp_perps = pysurvey_distance(r0,other_gals)

            # Calc para and perp ``our" way. 
            elif args.lado==False and args.oned==False and args.pysurvey==False:
                temp_paras,temp_perps = our_para_perp(r0,other_gals)

            # Calc Lado's way
            elif args.lado==True and args.oned==False and args.pysurvey==False:
                temp_paras,temp_perps = lado_para_perp(r0,other_gals)

            # Calc just the 1D
            elif args.lado==False and args.oned==True and args.pysurvey==False:
                temp_paras,temp_perps = one_dimension(r0,other_gals)

            paras[indexlo:indexhi] = temp_paras
            perps[indexlo:indexhi] = temp_perps

            indexlo = indexhi
        
        # Histogram the values.
        hist=plt.hist2d(perps[0:indexhi],paras[0:indexhi],bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
        tot_freq += hist[0]
        
        # Mirror the negative perps
        #hist=plt.hist2d(-1*perps[0:indexhi],paras[0:indexhi],bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
        tot_freq += hist[0]

        indexlo=0
        indexhi=0

        del hist

        print tot_freq.sum()
   
    print 'Point:'
    print tot_freq[199,101]
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
    np.savetxt(outfilename,tot_freq)



if __name__=='__main__':
    main()

# Two-Point Correlation Function
"""
DD=np.loadtxt('DD_1D_tot.dat',dtype='float')
DR=np.loadtxt('DR_1D_tot.dat',dtype='float')
RR=np.loadtxt('RR_1D_tot.dat',dtype='float')


DD = DD.transpose()
RR = RR.transpose()
DR = DR.transpose()

#DD+=np.flipud(DD)
#DR+=np.flipud(DR)
#RR+=np.flipud(RR)

#DDnew1=np.rot90(DD)
#DRnew1=np.rot90(DR)
#RRnew1=np.rot90(RR)
#DDnew2=np.rot90(DDnew1)
#DRnew2=np.rot90(DRnew1)
#RRnew2=np.rot90(RRnew1)

#DD+=DDnew2
#DR+=DRnew2
#RR+=RRnew2


'''
########## Bin Reduction ##########
def binred(array):
    import numpy as np
    binsnow=200
    binsneeded=100
    global newbins1
    global newbins
    newbins1=[]
    newbins=[]
    for i in range(0,200,2):
        #Adding columns
        global val1
        val1=(array[:,i]+array[:,i+1])
        newbins1.append(val1)
    newbins1=np.array(newbins1)
    return newbins1
    for j in range(0,200,2):
        global val
        val=(newbins1[j,:]+newbins1[j+1,:])
        newbins.append(val)
        return newbins
    newbins=np.array(newbins)    
                    
binred(DD)      
#print newbins        
#print val    
''' 




ndata=400000
nrand=1000000

#print DD.shape

#print sum(sum(DD))
#print sum(sum(DR))
#print sum(sum(RR))
'''
# Rebin
DDnew = np.zeros((100,100))
DRnew = np.zeros((100,100))
RRnew = np.zeros((100,100))


for i in range(0,200,2):
    for j in range(0,200,2):
        DDnew[i/2][j/2] = DD[i][j] + DD[i+1][j] + DD[i][j+1] + DD[i+1][j+1]
        DRnew[i/2][j/2] = DR[i][j] + DR[i+1][j] + DR[i][j+1] + DR[i+1][j+1]
        RRnew[i/2][j/2] = RR[i][j] + RR[i+1][j] + RR[i][j+1] + RR[i+1][j+1]
DD = DDnew
DR = DRnew
RR = RRnew

'''

DD /=(ndata**2-ndata)/2.
DR /=(nrand*ndata)/1.
RR /=(nrand**2-nrand)/2.
theta = (DD - 2*DR + RR)/RR

#R^2 WEIGHTING

nbins=400
rangeval=300

# Correct for little h
rangeval *= 0.7

#R Values
'''
for i in range(nbins):
    for j in range(nbins):
        r2=((nbins/2)-i)**2 + (j-(nbins/2))**2
        theta[i][j] *= r2
'''

plt.figure(figsize=(8,8))


#extent=
#plot=plt.imshow(theta)
extent=[-rangeval,rangeval,-rangeval,rangeval]

plt.subplot(2,2,1)
a=plt.imshow(DD,extent=extent)
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title('DD')

plt.subplot(2,2,2)
b=plt.imshow(RR,extent=extent)
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title('RR')

plt.subplot(2,2,3)
c=plt.imshow(DR,extent=extent)
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title('DR')

### Mirror Over the X-Axis #### 
#nbins = 100
newtheta= np.zeros((nbins,nbins))
newtheta += theta
for i in range(0,nbins):
    newtheta[i] += theta[(nbins-1)-i]

plt.subplot(2,2,4)
d=plt.imshow(newtheta,extent=extent,norm=mpl.colors.LogNorm(vmin=0.03,vmax=0.2))
plt.colorbar(d)
plt.xlabel(r'$r_\perp (h^{-1}$Mpc)')
plt.ylabel(r'$r_\parallel (h^{-1}$Mpc)')
plt.title(r'$\xi$')

#plt.ylim(0,1)
#plt.xlim(0,30)
#plt.xlabel('Distance (Mpc)')
#plt.ylabel('Theta')
#plt.title('DR10 Correlation Estimator with 20000 Galaxies')

plt.tight_layout()

plt.show()

if 1:
    w1D = theta.transpose()[200]
    w1D2= w1D[200:]
    xvals = np.linspace(-rangeval,rangeval,nbins/2)
    print w1D,xvals
    plt.figure()
    plt.plot(xvals,w1D2,'o')
    plt.xlim(0,200)
    plt.ylim(0,500)
    plt.xlabel('r')
    plt.ylabel('r^2 Xi')
    plt.title('1D Correlation Function (With r^2)')
plt.show()

#np.savetxt('ourxi.dat',newtheta)
"""


