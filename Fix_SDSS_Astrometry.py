#!/usr/bin/env python
import numpy as np
import scipy as sp
import sys, os, os.path, time, gc, glob
from astropy.table import Table
from astropy.io import ascii
from astropy.coordinates import SkyCoord
import astropy.io.fits as fits
import astropy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.wcs import WCS

# Things we need to manually input
date      = '20140321'                                              # Night we are doing
files     = [60,61,62,63,86,87,88,89,99,100,101,102]           # List of SDSS calibration files numbers
ra, dec  = '08h34m42.3s', '-03d31m34.7s'                            # RA and DEC of calibraion field (this reduces how long Astrometry.net takes to solve)

# This is to convert coords to degrees later
c = SkyCoord(ra, dec)

# Things we don't need to change
top         = '/mnt/Resources/perseus/CTIO_Data/' # path to CTIO data
reducedpath = top+date+'/Reduced_Data/'           # path to reduced CTIO data

######################################################################################################### Run Astrometry.net
# Run the command
for File in files:
    infile = reducedpath+date+'.{:0>3}.reduced.fits'.format(File)
    print(infile)
    command = "solve-field --no-plots --no-verify --ra %s --dec %s --radius 1 %s"%(c.ra.degree, c.dec.degree, infile)
    os.system(command)
#########################################################################################################
