#!/usr/bin/env python

"""Generate the CUTE input files from the SDSS catalogs and randoms posted at
http://data.sdss3.org/sas/dr11/boss/lss/

"""

from __future__ import division, print_function

import os
import sys
import argparse
import logging as log

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--dr', type=str, default='dr11', help='Specify the SDSS data release.')
    parser.add_argument('--parse', action='store_true', help='Parse the input datafiles.')
    parser.add_argument('--docute', action='store_true', help='Run CUTE.')
    parser.add_argument('--qaplots', action='store_true', help='Generate QAplots.')

    args = parser.parse_args()

    # Set up the logger; basic error checking.
    #log = logging.getLogger()

    key = 'LSS_BOSS'
    if key not in os.environ:
        log.fatal('Required ${} environment variable not set'.format(key))
        return 0

    drdir = os.path.join(os.getenv('LSS_BOSS'), args.dr)
    datafile = os.path.join(drdir, args.dr+'_cmass.dat')
    randomfile = os.path.join(drdir, args.dr+'_cmass_random.dat')
    outfile = os.path.join(drdir, 'dr11_2pt_rad.dat')
    paramfile = os.path.join(drdir, 'dr11_rad.param')

    # Parse the input data and write out CUTE-compatible files.
    if args.parse:
        allspecz = fits.getdata(os.path.join(drdir, 'galaxy_DR11v1_CMASS_North.fits.gz'), 1)
        keep = np.where((allspecz['Z']>0.43)*(allspecz['Z']<0.7))[0]

        specz = allspecz[keep]
        ngal = len(keep)

        data = np.zeros((ngal,4))
      
	data[:,0] = specz['RA']
	data[:,1] = specz['DEC']
	data[:,2] = specz['Z']
	data[:,3] = specz['WEIGHT_FKP']*specz['WEIGHT_SYSTOT']*(specz['WEIGHT_NOZ']+specz['WEIGHT_CP']-1)
        log.info('Writing {}'.format(datafile))
	np.savetxt(datafile, data)
	
	ra, dec, z, ipoly, wboss, wcp, wzf, veto = \
	  np.loadtxt(os.path.join(drdir, 'mock_random_DR11_CMASS_N_PTHALOS_ir4001.dat'), unpack=True)
	keep = np.where(veto==1)[0]
	nobj = len(keep)
	
	rand = np.zeros((nobj,4))
	rand[:,0] = ra[keep]
	rand[:,1] = dec[keep]
	rand[:,2] = z[keep]
	rand[:,3] = wcp[keep]+wzf[keep]-1
	
        plt.figure()
        plt.plot(data[:,0],data[:,1],'bo')
        plt.show()

        log.info('Writing {}'.format(randomfile))
	np.savetxt(randomfile, rand)

    if args.docute:
        # Do stuff; write paramfile; call cute using os.system()
        pfile = open(paramfile,'w')
        
        pfile.write('data_filename= '+datafile+'\n')
        pfile.write('random_filename= '+randomfile+'\n')
        pfile.write('input_format= 2\n')
        pfile.write('mask_filename= junk\n')
        pfile.write('z_dist_filename= junk\n')
        pfile.write('output_filename= '+outfile+'\n')
        pfile.write('num_lines= all\n')
        pfile.write('corr_type= radial\n')
        pfile.write('corr_estimator= LS\n')
        pfile.write('np_rand_fact= 8\n')
        pfile.write('omega_M= 0.3\n')
        pfile.write('omega_L= 0.7\n')
        pfile.write('w= -1\n')
        pfile.write('log_bin= 0\n')
        pfile.write('n_logint= 10\n')
        pfile.write('dim1_max= 75.0\n')
        pfile.write('dim1_nbin= 75\n')
        pfile.write('dim2_max= 75.0\n')
        pfile.write('dim2_nbin= 75\n')
        pfile.write('dim3_min= 0.4\n')
        pfile.write('dim3_max= 0.7\n')
        pfile.write('dim3_nbin= 1\n')
        pfile.write('radial_aperture= 1\n')
        pfile.write('use_pm= 1\n')
        pfile.write('n_pix_sph= 2048\n')
     
        pfile.close()
        #os.system('CUTE-noweights '+paramfile)
        os.system('CUTE '+paramfile)
    
    if args.qaplots:
        # Make rockin' plots and write out.
        cutedata=np.loadtxt(outfile)
        for ii in range(0,len(cutedata[0])): # data -> cutedata
           globals()['matrix{0}'.format(ii)]=cutedata[:,ii]
        plt.figure()
        plt.loglog(matrix0,matrix1,'bo')
        plt.show()
        pass
           
if __name__ == "__main__":
    main()
