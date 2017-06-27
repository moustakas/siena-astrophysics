#!/usr/bin/env python

# J. Moustakas - 2014 Jun 11
#
# Convert the binary snapshot files written out by Gadget2 to ASCI
# text files which can be read by blendergadget.py

import os
import glob
import pynbody as pyn
import numpy as np
from optparse import OptionParser
from scipy import ndimage

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-d", "--directory", action="store", dest="directory", 
                      help="Directory of snapshot files", default="./")
    (opts, args) = parser.parse_args()

    print('Selected snapshot directory is '+opts.directory)

    if len(args)==0:
        snapfiles = sorted(glob.glob(opts.directory+'snapshot_*'))
    else:
        snapfiles = args
    nsnap = len(snapfiles)
    print('Processing '+str(nsnap)+' snapshot files.')
#   print(snapfiles)
    
# make the output directory
    txtdir = opts.directory+'/txt'
    if not os.path.exists(txtdir):
        print('Making output directory '+txtdir)
        os.makedirs(txtdir)

    for sfile in snapfiles:
        print('Reading '+sfile)

        sim = pyn.load(sfile)
        xyz_dm = sim.dm['pos']
        xyz_stars = sim.stars['pos']
        xyz_gas = sim.gas['pos']

# shift to the center-of-mass of the stars
        mass_stars = sim.stars['mass']
#       mass_gas = sim.gas['mass']

        com_stars = np.average(xyz_stars,axis=0,weights=mass_stars)
#       com = ndimage.measurements.center_of_mass(xyz_stars)
        xyz_stars_new = xyz_stars-com_stars
        xyz_gas_new = xyz_gas-com_stars
        
# write the positions out        
        starfile = txtdir+'/'+os.path.basename(sfile)+'_stars.txt'
        gasfile = txtdir+'/'+os.path.basename(sfile)+'_gas.txt'
        print('Writing '+starfile+', '+gasfile)

        np.savetxt(starfile,xyz_stars_new,fmt='%.3e')
        np.savetxt(gasfile,xyz_gas_new,fmt='%.3e')
