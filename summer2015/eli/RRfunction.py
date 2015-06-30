@profile
def RR():
    import numpy as np
    import math
    import scipy.spatial
    from astropy.cosmology import FlatLambdaCDM
    import astropy.io
    import matplotlib.pylab as plt
    import time
    t0=time.time()
    r=np.loadtxt('cmass_dr10_north_randoms_ir4500.v10.1.release.txt')
    ra=r[:,0]
    dec=r[:,1]
    z=r[:,2]


    ##### Z Cut #####

    zcut=z>0.43
    zcutnew=z[zcut]

    zcut1=z<0.7
    zcutnew1=z[zcut1]
    tot=zcut*zcut1
    totrand=r[tot]
    del tot
    del r
    ########### Distances ################

    nrands_for_calculation = 250000

    np.random.seed(1)

    a=np.arange(0,len(totrand))
    np.random.shuffle(a)
    sample=totrand[a[0:nrands_for_calculation]]
    del a
    del totrand
    ###### Distances (Para and Perp) ########

    print 'Conversions'
    # Comoving Distances
    cosmo=FlatLambdaCDM(H0=70,Om0=0.3)
    comdist=cosmo.comoving_distance(sample[:,2])

    # Converting RA and Dec to Radians

    RArad=(sample[:,0])*((math.pi)/180)
    Decrad=((math.pi)/2)-((sample[:,1])*((math.pi)/180))

    # Converting to Cartesian Coordinates

    x=comdist*np.sin(Decrad)*np.cos(RArad)
    y=comdist*np.sin(Decrad)*np.sin(RArad)
    z=comdist*np.cos(Decrad)
    coordsa=np.column_stack((x,y,z))
    del RArad
    del Decrad
    del x
    del y
    del z
    del comdist
    del sample
    # Distances

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
    nchunks = nrands_for_calculation/chunk_size
    nbins=200
    rangeval=300


    tot_freq = np.zeros((nbins,nbins)) 

    for j in xrange(nchunks):
        lo = j*chunk_size
        hi = (j+1)*chunk_size
        print "Performing calculations for RR %d chunk: %d-%d" % (j,lo,hi)

        paras = []
        perps = []

        for i in range(lo,hi):
            if i!=ngals-1:
                # First compute R_LOS and dR
               
                r1=coordsa[i]
                
                R_LOS1 = (r1 + coordsa[i+1:])/2.
                
                dR1 = coordsa[i+1:] - r1
                
                R_LOS_mag1 = mag(R_LOS1)

                # Dot product
                
                R_para1 = (dR1[:,0]*R_LOS1[:,0] + dR1[:,1]*R_LOS1[:,1] + dR1[:,2]*R_LOS1[:,2])/R_LOS_mag1
                
                dR_mag1 = mag(dR1)
                # Make use of the Pythagorean theorem
                
                R_perp1 = np.sqrt(dR_mag1*dR_mag1 - R_para1*R_para1)
               
                
                paras += R_para1.tolist()
                
                perps += R_perp1.tolist()
                

        
        hist=plt.hist2d(perps,paras,bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
        
        tot_freq += hist[0]
        
        # Mirror the negative perps
        hist=plt.hist2d(-1*np.array(perps),paras,bins=nbins,range=((-rangeval,rangeval),(-rangeval,rangeval)))
        tot_freq += hist[0]
        
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
    np.savetxt('RRtest2d2.txt',tot_freq)
    t1=time.time()
    tottime=t1-t0
    totmin=tottime/60
    tothr=totmin/60
    print 'This code took %f hours to run' %tothr
RR()
