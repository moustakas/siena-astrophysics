#!/usr/bin/python

from __future__ import division, print_function

import os
import sys
import argparse
import glob
import logging as log
import numpy as np

# create directory and environment variable
CosmoMC = os.path.join(os.getenv('COSMOMC'))
outdir = os.path.join(Cosmo, 'cosmoMCout')
paramfile = os.path.join(Cosmo, 'param')

