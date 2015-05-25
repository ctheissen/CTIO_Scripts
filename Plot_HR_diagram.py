#!/usr/bin/env python
import numpy as np
import sys, os, os.path, time, gc, glob
from astropy.table import Table
from astropy.io import ascii
import photutils as ph
import astropy.io.fits as fits
import astropy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.wcs import WCS

# Things we need to input
date        = '20140128'  # Night of the observations

# Things we don't need to input
top         = '/mnt/Resources/perseus/CTIO_Data/'
reducedpath = top+date+'/Reduced_Data/'

# Plotting params
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
plt.rc('text', usetex=True)
plt.rc('legend', fontsize=10)
plt.rc('axes', labelsize=14)

######################################################################################################### Read in the file
#print(os.getcwd())
#print(glob.glob(os.getcwd()+'/*.csv'))
#sys.exit()
data0 = np.genfromtxt(glob.glob(os.getcwd()+'/*.csv')[0], delimiter=',', usecols=(0,1,9,10,15,16,17,18,20), missing_values = -9999, filling_values = -9999,
                      dtype=[('Obj', np.str_, 20),('Num', np.int),('zmagerr', '<f8'), ('pi', '<f8'), ('Mr', '<f8'), ('Mr_err', '<f8'),
                             ('rz', '<f8'), ('rz_err', '<f8'), ('sub', '<i8')])
print(data0.dtype.names)
data1 = data0[np.where( (data0['Num']==1) & (data0['sub']!=1) & (data0['zmagerr']!=-9999) & (data0['pi']!=-9999) )] # Remove bad data/binaries/subdwarfs
data2 = data0[np.where( (data0['Num']>1) & (data0['sub']!=1) & (data0['zmagerr']!=-9999) & (data0['pi']!=-9999) )] # Get just the binaries
data3 = data0[np.where( (data0['Num']==1) & (data0['sub']==1) & (data0['zmagerr']!=-9999) & (data0['pi']!=-9999) )] # Get the subdwarfs
print(len(data3), len(data2), len(data1))

rz   = data1['rz']
r    = data1['Mr']
rzUn = data1['rz_err']
rUn  = data1['Mr_err']

Bochanski = lambda x: 5.190 + 2.474*x + 0.4340*x**2 - 0.08635*x**3
West = lambda x: 5.190 + 2.474*x + 0.4340*x**2 - 0.08635*x**3
xp = np.linspace(rz.min(),rz.max())

plt.figure(101, figsize=(8,5))
#scat1, = plt.errorbar(rz, r, xerr = rzUn, yerr = rUn, ls='None', c='b')
plt.errorbar(data1['rz'], data1['Mr'], xerr = data1['rz_err'], yerr = data1['Mr_err'], ls='None', c='b', fmt='o', markersize=2, markeredgecolor='b')
plt.errorbar(data2['rz'], data2['Mr'], xerr = data2['rz_err'], yerr = data2['Mr_err'], ls='None', c='m', fmt='o', markersize=2, markeredgecolor='m')
plt.errorbar(data3['rz'], data3['Mr'], xerr = data3['rz_err'], yerr = data3['Mr_err'], ls='None', c='r', fmt='o', markersize=2, markeredgecolor='r')
plt.plot(xp, Bochanski(xp), 'r--')
ylower, yupper = plt.ylim()
plt.ylim(yupper, ylower)
plt.xlabel(r'$r-z$')
plt.ylabel(r'$M_r$')
plt.legend([r'${\rm Clean}$', r'${\rm Binaries}$', r'${\rm Subdwarfs}$', r'${\rm Bochanski}$'], numpoints=1)
plt.savefig('Plots/HR.png', dpi=600)
plt.savefig('Plots/HR.pdf')
plt.show()
