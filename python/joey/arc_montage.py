#! /usr/bin/env python

import glob
import numpy as np

def arc_montage():
    path = '/home/obsastro1/clash/'
    all_arcs = glob.glob(path + 'macs1206-arc*-im*.png')
    narcs = len(all_arcs)

    gal = np.zeros(narcs)
    for ii in range(narcs):
        gal[ii] = float(all_arcs[ii].split('-')[1][3:])
        print gal[ii]

if __name__ == '__main__':
    arc_montage()
