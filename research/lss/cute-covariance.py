#!/usr/bin/python

from __future__ import division, print_function

import os
import sys
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
import pylab

def main():
    
    # argument parsing
    parser = argparse.ArgumentParser()

    parser.add_argument('--dr', type=str, default='dr11', help='Specify the SDSS data release.')
    parser.add_argument('--type', type=str, default='3D_rm', help='Specify the correlation type used.')
    parser.add_argument('--cov', action='store_true', help='Compute the covariace matrix.')

    args = parser.parse_args()

    # convenience variables
    datadir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'cuteout', args.type)
    allfiles = glob.glob(os.path.join(datadir, '*.dat'))
    outdir = os.path.join(os.getenv('LSS_BOSS'), args.dr, 'covariance')

    rad = []
    mu = []
    xi = []
    covariance = []
    ximeans = []
    vectors = []
    
    if args.cov:

        for ii in range(len(allfiles)):
            mu.append(np.loadtxt(allfiles[ii])[:,0])
            rad.append(np.loadtxt(allfiles[ii])[:,1])
            xi.append(np.loadtxt(allfiles[ii])[:,2])
            print(ii)

        xi = np.asarray(xi)
        radbar = np.mean(rad)
        mubar = np.mean(mu)

        for ii in range(len(xi)):
            ximeans.append(np.mean(xi[ii]))

        #mean = np.mean(ximeans)
        
        for ii in range(xi.shape[0]):
            if ii >= 1:
                solution = np.cov(xi[ii-1], xi[ii])
                matrix =  np.matrix(solution)
                covariance.append(matrix)

        for ii in range(len(covariance)):
            vector, value = np.linalg.eig(covariance[ii])
            vectors.append(vector)

        for ii in range(len(vectors)):
            plt.scatter(vectors[ii][0], vectors[ii][1])

        plt.show()
        pylab.pcolor(np.corrcoef(mu));pylab.colorbar();plt.show()
            
            
if __name__ == "__main__":
    main()
    
