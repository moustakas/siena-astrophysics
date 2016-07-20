===========================
Siena Large Scale Structure
===========================

Introduction
------------

This directory contains code developed by students and faculty at Siena College
to investigate large-scale structure and the baryonic acoustic oscillations
(BAO) feature using the 2-point correlation function.

The only dataset currently supported is the SDSS-III/BOSS Data Release 11 (DR11)
CMASS/North spectroscopic redshift catalog, although additional datasets
(including the BOSS/DR12 and simulated DESI catalogs) will be added in the
future.

The code we use to compute the correlation function is the public `CUTE`_ code
(see also the `CUTE Github repository`_), which is written in C.  For the random
catalogs we adopt the `PTHalo mocks`_ from Manera et al. 2013, 2015.

.. _`CUTE`: http://members.ift.uam-csic.es/dmonge/CUTE.html

.. _`CUTE Github repository`: https://github.com/damonge/CUTE

.. _`PTHalo mocks`: http://www.marcmanera.net/mocks


Install Required Dependencies
-----------------------------

First, clone the ``siena-astrophysics`` repository somewhere convenient (or just
grab the script ``cute-2ptcor.py``, as there are no other internal dependencies)
and, for convenience, create an alias (which you may choose to put in your
startup file):

.. code:: bash
          
          git clone git@github.com:moustakas/siena-astrophysics.git
          alias cute-2ptcor=`python /path/to/sienarepo/research/lss/cute-2ptcor.py'

Next, download `CUTE`_, read or skim the documentation, compile it, and put the
executable in your path.  When compiling be sure to compile *with* the
``-D_WITH_WEIGHTS`` option.  Note:

- ``CUTE`` can be compiled either with or without MPI support, depending on the
machine you're planning to use to run it.

- We do not yet support the ``CUDA`` version of ``CUTE``, so that version does
  not have to be compiled.

Finally, the code also depends on various standard Python libraries like
``numpy`` and ``matplotlib``.

Download Required Files
-----------------------

First, create a top-level directory where you want to keep all the input/output
files and then create an environment variable called ``$LSS_CUTE`` pointing to
that directory.  For example, in bash,

.. code:: bash

        mkdir /path/to/lss-cute
        export LSS_CUTE=/path/to/lss-cute

Next, create all the sample-specific subdirectory we need.  Here, we give
instructions for handling just the CMASS/North DR11 sample.

.. code::

        mkdir $LSS_CUTE/dr11_cmass_north

Third, download the binary FITS table containing the BOSS/DR11 large-scale
structure catalog (approximately 100M):

.. code::
          
        cd $LSS_CUTE/dr11_cmass_north
        wget http://data.sdss3.org/sas/dr11/boss/lss/galaxy_DR11v1_CMASS_North.fits.gz

Finally, download and unpack the PTHalo mock catalogs (approximately 4.3G):

.. code:: bash
          
        mkdir $LSS_CUTE/dr11_cmass_north/randoms
        cd $LSS_CUTE/dr11_cmass_north/randoms
        wget http://data.sdss3.org/sas/dr11/boss/lss/dr11_pthalos_mocks/mock_random_DR11_CMASS_N_PTHALOS_allmocks.tar.gz
        tar xzvf mock_random_DR11_CMASS_N_PTHALOS_allmocks.tar.gz

At this point you may delete the
``mock_random_DR11_CMASS_N_PTHALOS_allmocks.tar.gz`` file to save space,
otherwise you're ready to go!

Examples
--------

``cute-2ptcor`` is a simple command-line script which interfaces with CUTE in a
convenient way (for example, by generating the required parameter files
on-the-fly).  It provides a non-exhaustive set of options for computing the
correlation function, and more options can be added as needed.  To get the full
set of available options 

.. code:: bash

          usage: cute-2ptcor.py [-h] [--dr DR] [--parse] [--docute] [--qaplots]
                      [--corrtype CORRTYPE] [--cosmo COSMO]

          optional arguments:
            -h, --help           show this help message and exit
            --dr DR              Specify the SDSS data release.
            --parse              Parse the input datafiles (do just once).
            --docute             Run CUTE.
            --qaplots            Generate QAplots.
            --corrtype CORRTYPE  Specify correlation type (monopole|3D_ps|3D_rm).
            --cosmo COSMO        Adopted cosmology (1|2)


Compute and Plot the Monopole
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first time the code is run it reads in the spectroscopic catalog as well as
the random catalogs (or a subset thereof) and parses them.

.. code:: python
          stuff

You can also choose from among two possible cosmologies:


Compute the 2D Correlation Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, compute the 2D correlation function in pi-mu space.

.. code:: python
          stuff

Compare with the literature.




Contributors
------------

- Elijah Beaudin (2016)
  
- Kevin Napier (2018)
  
- Prof. John Moustakas (Physics)
  
- Prof. Matt Bellis (Physics)
