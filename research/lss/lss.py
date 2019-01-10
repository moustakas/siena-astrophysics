#! /usr/bin/env python

import os, pdb
import numpy as np
from mpi4py import MPI
from nbodykit.lab import *
from nbodykit import setup_logging

cosmo = cosmology.Cosmology(h=0.7).match(Omega0_m=0.31)
datadir = os.path.join(os.getenv('IM_DATA_DIR'), 'sdss', 'dr12')

setup_logging()

def subsample2_data(randoms=False, nsample=1, clobber=False):
    """Read and subsample the data or randoms LSS catalog, for speed.
    
    Here, we select objects with z=0.4-0.6 (see Alam+17, Table 2).
    
    """
    import fitsio

    if randoms:
        dataset = 'DR12v5-randoms'
        columns = ['RA', 'DEC', 'Z', 'WEIGHT_FKP',]
        infile = os.path.join(datadir, 'random0_DR12v5_CMASSLOWZTOT_North.fits.gz')
        outfile = os.path.join(datadir, 'random0_DR12v5_CMASSLOWZTOT_North')
    else:
        dataset = 'DR12v5'
        columns = ['RA', 'DEC', 'Z', 'WEIGHT_FKP', 'WEIGHT_CP', 'WEIGHT_NOZ', 'WEIGHT_SYSTOT',]
        infile = os.path.join(datadir, 'galaxy_DR12v5_CMASSLOWZTOT_North.fits.gz')
        outfile = os.path.join(datadir, 'galaxy_DR12v5_CMASSLOWZTOT_North')
    try:
        data = BigFileCatalog(outfile, dataset='samples')
    except:
        clobber = True
        
    if clobber: # convert to bigfile if not yet.
        # MPI stuff here
        data = ArrayCatalog(fitsio.read(infile, columns=columns))
        # if data doesn't fit into memory, use FITSCatalog(infile)
        data.save(outfile, columns=columns, dataset='samples')
    
    data = BigFileCatalog(outfile, dataset='samples')#, comm=None)
    
    if nsample is not None:
        data = data[::nsample]
        
    keep = (data['Z'] > 0.4) & (data['Z'] < 0.6)
    print(keep.sum().compute())
    out = data[keep]
        
    out['Position'] = transform.SkyToCartesian(out['RA'], out['DEC'], out['Z'], cosmo=cosmo)
    
    return out

def get_nofz(totarea=9329, north=True):
    """Build n(z) for the data.  The total area is from Alam+17, Table 2.
    
    """
    from scipy.interpolate import InterpolatedUnivariateSpline
    
    if north:
        area = totarea * 5.3 / 7.3
    else:
        area = totarea * 2.0 / 7.3
    
    # compute n(z) from the randoms
    fsky = area / ( 4 * np.pi * (180 / np.pi)**2 ) # sky fraction from Alam+17
    zhist = RedshiftHistogram(random, fsky, cosmo, redshift='Z')
    
    # re-normalize to the total size of the data catalog
    alpha = 1.0 * data.csize / random.csize
    print('Renormalization factor = {:.5f}'.format(alpha))
    zhist.nbar *= alpha
    
    # Compute the interpolator we'll use below.
    nofz = InterpolatedUnivariateSpline(zhist.bin_centers, zhist.nbar)
    
    return zhist, nofz

#data['NZ'] = nofz(data['Z'])
#random['NZ'] = nofz(random['Z'])
#data['FKPWeight'] = 1.0 / (1 + data['NZ'] * 2e4)
#random['FKPWeight'] = 1.0 / (1 + random['NZ'] * 2e4)

#fkp = FKPCatalog(data, random)
#mesh = fkp.to_mesh(Nmesh=256, nbar='NZ', comp_weight='Weight', fkp_weight='FKPWeight')

#zhist, nofz = get_nofz()



data = subsample2_data(nsample=1, clobber=False)
random = subsample2_data(nsample=1, randoms=True, clobber=False)

data['Weight'] = (data['WEIGHT_SYSTOT'] * (data['WEIGHT_NOZ'] + data['WEIGHT_CP'] - 1)).compute()
random['Weight'] = np.ones(len(random))

pimax = 2
edges = np.array([0.5, 1, 2, 3])

xi = SurveyData2PCF('projected', data[::100], random[::100], edges, cosmo=cosmo, pimax=pimax,
                    redshift='Z')#, weight='Weight', show_progress=True)#, **{'nthreads': 4})
# xi.save()


#pdb.set_trace()
