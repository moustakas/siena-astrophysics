# An Exploration of the Stellar Mass Function of BCGs

* This directory is dedicated to the study of the masses of BCGs in galaxy 
clusters. 

* ***NOTE:*** My own (very basic) CLF calculator is under lumfunc.ipynb. The variable names suggest that the mass function is being studied in the notebook. However, this is just a consequence of my laziness when making modifications; relic of the notebook's development. It is in semi-functioning condition, returning the expected lognormal distribution for ```latex $\Phi (L)$``` versus ```latex $L (L_{\odot})$```. 

* Reddick's code and files are available on NYX at: `/global/work/projects/redmapper/csmf`

* <a href="https://www.evernote.com/shard/s730/sh/982bd8b9-3c15-4fe1-be5c-5068f7995aa3/936bc9c36eb37628" target="_blank">Here are Reddick's notes on the code.</a>


##### A few notes on the use of the code:

* The codes are all written in Python 2. 
* `chtolib` ***needs*** to be in the `redmapper` directory at the moment.
* The code requires `pyfits` (installed using pip)
* The code requires `esutil` (installed using pip)
* The code requires `healpy` (installed using pip)


##### Calculating the Conditional Stellar Mass Function (CSMF):

* To calculate the CSMF, use the `redm_full.py` script in the `redmapper` directory. Use the following in the command line:`python redm_full.py paramfile`

* The paramfile must be an ascii file. It specifies input files, an output directory, and some optional flags. For an example, see `paramfile_dr8.dat` in the `redmapper` directory. The input files are asignes to rows of the following names:
  ** cluster_file (central galaxies catalog)
  ** member_file (satellites catalog)
  ** kcorr_file
  ** cindex_file (indexes centrals to the members list)

* ***NOTE:*** The code does **not** internally account for changes in effective area with redshift


##### Fitting the CSMF:

* The first step before doing anything else is to calculate the $\lambda (M)$ relationship.
  ** This is not yet an automated process. You will have to run some variation of the following bit of code in an `iPython notebook`.
```python
import numpy as np
import mass_matching
import pyfits
cat = pyfits.open("/nfs/slac/g/ki/ki19/des/erykoff/clusters/process/dr8_redmapper_v6.2/run_ubermem/dr8_run_redmapper_v6.2_ubermem_lgt5_catalog.fit") #This reads in the clusters catalog -- note members are not needed
cat = cat[1].dat
zmin = np.array([0.1, 0.15, 0.2, 0.25, 0.3]) #Min ranges for z-bins for matching
zmax = np.array([0.15, 0.2, 0.25, 0.3, 0.33]) #Max of the range -- so, the first bin is z=[0.1,0.15]
zlabel = ['0.125','0.175','0.225','0.275','0.315'] #Text labels -- used to get midpoint mass functions
area_bin = np.repeat(10405.,5) #Effective area of each redshift bin; is the same for most, but included for, e.g., sva, where area is not constant
outdir = "/nfs/slac/g/ki/ki10/rmredd/redmapper_data/dr8_zlambda_v6.2_ubermem/" #Output directory for the lambda(mass) results
param, param_err = mass_matching.run_matching_zbin_and_print(cat,zmin,zmax,zlabel,area_bin,outdir)
```

* You may now runt he function `run_zev_fit.py`: `python run_zev_fit.py A_lm_z.dat constant 0.1 0.33 0.9 10` along with the output files, to obtain the redshift evolution parameters. You must perform this twice. The second time, you should instead pass the argument, `lnlm0_z.dat power_log 0.1 0.33 3. 0.85 10`.
