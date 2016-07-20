=====================
Large Scale Structure
=====================

Introduction
------------

This directory contains code to investigate large-scale structure and the BAO
feature using the 2-point correlation function.  The dataset we currently
support is the SDSS-III/BOSS Data Release 11 CMASS-North spectroscopic redshift
catalog, and the driver code we use to compute the correlation function is the
public CUTE code, which is written in C++.  For the random catalogs we adopt the
PTHalo mocks from Manera et al. 2013, 2015.

Install Required Dependencies
-----------------------------

First, clone the '''siena-astrophysics''' repository (or just grab
'''research/lss/cute-2ptcor.py''', as there are no other internal dependencies).

Download CUTE, skim the documentation, compile it, and put the executable in
your path.  When compiling be sure to compile *with* the -D_WITH_WEIGHTS option.
Note: we do not use the CUDA version of the code and CUTE can be compiled either
with or without MPI support.

The code also depends on various standard Python libraries like '''numpy''' and
'''matplotlib'''.

Download Required Files
-----------------------

First, create a top-level directory where you want to keep all the input/output
files and then create an environment variable called LSS_BOSS pointing to that
directory.  For example, in bash,

'''
export LSS_BOSS=/path/to/lss/dr11
'''

Next, download the binary FITS table containing the BOSS/DR11 large-scale
structure catalog:
'''
cd $LSS_BOSS
wget http://data.sdss3.org/sas/dr11/boss/lss/galaxy_DR11v1_CMASS_North.fits.gz
'''

Finally, download and unpack the PTHalo mock catalogs:

'''
mkdir $LSS_BOSS/randoms
cd $LSS_BOSS/randoms
wget http://data.sdss3.org/sas/dr11/boss/lss/dr11_pthalos_mocks/mock_random_DR11_CMASS_N_PTHALOS_allmocks.tar.gz
tar xzvf mock_random_DR11_CMASS_N_PTHALOS_allmocks.tar.gz
'''

You may delete the '''.tar.gz''' file to save space if you'd like, otherwise
you're ready to go!

Compute and Plot the Monopole
-----------------------------

To compute the monopole 
