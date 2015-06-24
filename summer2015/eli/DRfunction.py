@profile
def DR():
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
    t0=time.time()

    # Opening FITS file (SDSS Data) 'a' Argument = dr10cmassnorth.fits
    print 'Reading in FITS Data'
    infilename1 = sys.argv[1]
    hdulist1=fits.open(infilename1)
    hdulist1.info()
    h1=hdulist1[1]
    data1=h1.data

    # Opening txt file (Mocks) 'b'
    print 'Reading in Text File'
    r=np.loadtxt('cmass_dr10_north_randoms_ir4500.v10.1.release.txt')
    z=r[:,2]

    print "Read in text file......"

    print "Making cuts......"
    # Txt Z Cut
    zcut=z>0.43
    zcutnew=z[zcut]

    zcut1=z<0.7
    zcutnew1=z[zcut1]

    tot=zcut*zcut1
    totrand=r[tot]
    del tot
    del r

    # Randomizing a Sample of SDSS Data
    ngals_for_calculation = 200000
    nrands=250000
    np.random.seed(1)

    a=np.arange(0,len(data1))
    np.random.shuffle(a)
    samplea=data1[a[0:ngals_for_calculation]]
    del h1
    del data1
    del a
    # Randomizing a sample of Mock Data
    b=np.arange(0,len(totrand))
    np.random.shuffle(b)
    sampleb=totrand[b[0:nrands]]
    del b
    del totrand
    ##### Finding Values for Spherical Coordinates ####

    # Comoving Distances
    cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
    comdista=cosmo.comoving_distance(samplea['Z'])
    comdistb=cosmo.comoving_distance(sampleb[:,2])

    ##### Converting RA and Dec to Radians #####

    # SDSS
    RArada=(samplea['PLUG_RA'])*((math.pi)/180)
    Decrada=((math.pi)/2)-((samplea['PLUG_DEC'])*((math.pi)/180))

    # Mock
    RAradb=(sampleb[:,0])*((math.pi)/180)
    Decradb=((math.pi)/2)-((sampleb[:,1])*((math.pi)/180))


    ##### Converting to Cartesian Coordinates #####

    # SDSS
    xa=comdista*np.sin(Decrada)*np.cos(RArada)
    ya=comdista*np.sin(Decrada)*np.sin(RArada)
    za=comdista*np.cos(Decrada)
    coordsa=np.column_stack((xa,ya,za))

    del RArada
    del Decrada
    del xa
    del ya
    del za
    del comdista
    del samplea
    # Mock
    xb=comdistb*np.sin(Decradb)*np.cos(RAradb)
    yb=comdistb*np.sin(Decradb)*np.sin(RAradb)
    zb=comdistb*np.cos(Decradb)
    coordsb=np.column_stack((xb,yb,zb))

    del RAradb
    del Decradb
    del xb
    del yb
    del zb
    del comdistb
    del sampleb
    print 'Finished with conversions! Now to calculate distances...'
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

    ngals = len(coordsa)
    chunk_size = 50
    nchunks = ngals_for_calculation/chunk_size
    nbins=200
    rangeval=300


    tot_freq = np.zeros((nbins,nbins)) 

    for j in xrange(nchunks):
        lo = j*chunk_size
        hi = (j+1)*chunk_size
        print "Performing calculations for DR %d chunk: %d-%d" % (j,lo,hi)

        paras = []
        perps = []
      
        for i,r1 in enumerate(coordsa[lo:hi]):
                # First compute R_LOS and dR
                R_LOS1 = (r1 + coordsb[:])/2
                dR1 = coordsb - r1
                R_LOS_mag1 = mag(R_LOS1)

                # Dot product
                
                R_para1 = (dR1[:,0]*R_LOS1[:,0] + dR1[:,1]*R_LOS1[:,1] + dR1[:,2]*R_LOS1[:,2])/R_LOS_mag1
                
                dR_mag1 = mag(dR1)
                # Make use of the Pythagorean theorem
                
                R_perp1 = np.sqrt(dR_mag1*dR_mag1 - R_para1*R_para1)
                #negR_perp1 = -1*R_perp1
                
                paras += R_para1.tolist()
                
                perps += R_perp1.tolist()
                #nperps1 += negR_perp1.tolist()
                #if i%(chunk_size/4)==0:
                    #print i

        #print len(paras)
        #print len(perps)
        #newperps1=np.concatenate((perps1,nperps1))
        #newparas1len=np.concatenate((paras1,paras1))

        #print 'Histogram1'

        
        hist=plt.hist2d(perps,paras,bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
        
        tot_freq += hist[0]
        
        # Mirror the negative perps
        hist=plt.hist2d(-1*np.array(perps),paras,bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
        tot_freq += hist[0]


        #print type(hist1[0])
        #frequ1=hist1[0]
        #plt.close()

        del paras
        del perps
        del hist

        print tot_freq
        
    #tot_freq[(nbins/2),(nbins/2)]=0
    print 'Final Plot'    
    #extent = [-rangeval,rangeval, -rangeval,rangeval]
    #fig = plt.figure()
    #axes = fig.add_subplot(1,1,1)
    print 'Imshow'
    #ret = axes.imshow(tot_freq,extent=extent,interpolation='nearest') #,origin=origin,cmap=cmap,axes=axes,aspect=aspect
    #plt.show()
    np.savetxt('DRtest2d1.txt',tot_freq)
    t1=time.time()
    tottime=t1-t0
    totmin=tottime/60
    tothr=totmin/60
    print 'This code took %f hours to run' %tothr
DR()
