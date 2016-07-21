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

- ``CUTE`` can be compiled either with or without MPI support, depending on the machine you're planning to use to run it.

- We do not yet support the ``CUDA`` version of ``CUTE``, so that version does not have to be compiled.

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

Next, create the sample-specific subdirectory we need, and download the binary
FITS table containing the BOSS/DR11 large-scale structure catalog (approximately
100M).  Here, we give instructions for handling just the CMASS/North DR11
sample.

.. code::

        mkdir $LSS_CUTE/dr11_cmass_north
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

``cute-2ptcor`` is a simple command-line script which interfaces with ``CUTE``
in a convenient way (for example, by generating the required parameter files
on-the-fly).  It provides a non-exhaustive set of options for computing the
correlation function and uses sensible defaults whenever possible; more options
can be added as needed.

The optional inputs to the code can be inspected by executing ``cute-2ptcor -h``
from the command line, or by calling the code without any inputs.  Here is the
current output:

.. code:: bash

	usage: cute-2ptcor.py [-h] [--sample SAMPLE] [--omegaM OMEGAM] [--w W]
	                      [--corrtype CORRTYPE] [--nrandom NRANDOM] [--docute]
	                      [--qaplots] [--clobber]
	
	optional arguments:
	  -h, --help           show this help message and exit
	  --sample SAMPLE      Dataset to use (currently only dr11_cmass_north is
	                       supported).
	  --omegaM OMEGAM      Omega_matter (note: Omega_Lambda = 1-Omega_Matter)
	  --w W                w parameter (choose w=-1.0 for cosmological constant)
	  --corrtype CORRTYPE  Specify correlation type (monopole|3D_ps|3D_rm).
	  --nrandom NRANDOM    Number of random catalogs to use (integer number|all)
	  --docute             Generate the individual correlation functions using
	                       CUTE.
	  --qaplots            Generate QAplots.
	  --clobber            Regenerate the parsed data/random files, even if they
	                       exist.


Compute and Plot the Monopole
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To compute and generate a plot of the simple monopole correlation function using
a single random catalog and all other defaults do:

.. code:: python

        cute-2ptcor --corrtype monopole --nrandom 1 --docute --qaplots

The first time the code is run it reads in the data and (needed) random catalogs
and parses them into the format read by ``CUTE``; in subsequent calls the code
only regenerate these catalogs if the optional input ``--clobber`` is set.

The output file is written to the ``qaplots`` directory.


Compute the 2D Correlation Function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Next, compute the 2D correlation function in pi-mu space.

.. code:: python


Compare with the literature...

Contributors
------------

- Elijah Beaudin (Class of 2016)
  
- Kevin Napier (Class of 2018)
  
- Prof. Matt Bellis (Siena College)

- Prof. John Moustakas (Siena College)
