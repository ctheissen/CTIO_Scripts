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
"""
data0 = np.genfromtxt(glob.glob(os.getcwd()+'/*.csv')[0], delimiter=',', usecols=(0,1,9,10,15,16,17,18,20,21), missing_values = -9999, filling_values = -9999,
                      dtype=[('Obj', np.str_, 20),('Num', np.int),('zmagerr', '<f8'), ('pi', '<f8'), ('Mr', '<f8'), ('Mr_err', '<f8'),
                             ('rz', '<f8'), ('rz_err', '<f8'), ('sub', '<i8')])
"""
data0 = np.genfromtxt(glob.glob(os.getcwd()+'/*.csv')[0], delimiter=',', usecols=(0,1,4,5,6,7,8,9,10,15,16,17,18,21,22,23), missing_values = -9999, filling_values = -9999,
                      dtype=[('Obj', np.str_, 20),('Num', np.int),('rmag', '<f8'),('rmagerr', '<f8'),('imag', '<f8'),('imagerr', '<f8'),('zmag', '<f8'),('zmagerr', '<f8'),
                             ('pi', '<f8'), ('Mr', '<f8'), ('Mr_err', '<f8'), ('rz', '<f8'), ('rz_err', '<f8'), ('sub', '<i8'), ('date', np.str_, 20),
                             ('person', np.str_, 20)])

print(data0.dtype.names)
data1 = data0[np.where( (data0['Num']==1) & (data0['sub']!=1) & (data0['rmag']!=-9999) & (data0['rmagerr']!=-9999) & (data0['pi']!=-9999) &
                        (data0['imag']!=-9999) & (data0['imagerr']!=-9999) & (data0['zmag']!=-9999) & (data0['zmagerr']!=-9999) )] # Remove bad data/binaries/subdwarfs
data2 = data0[np.where( (data0['Num']>1) & (data0['sub']!=1) & (data0['zmagerr']!=-9999) & (data0['pi']!=-9999) )] # Get just the binaries
data3 = data0[np.where( (data0['Num']==1) & (data0['sub']==1) & (data0['zmagerr']!=-9999) & (data0['pi']!=-9999) )] # Get the subdwarfs
print(len(data3), len(data2), len(data1))

################ First plot r-i vs. i-z
plt.figure(101, figsize=(20,10))
titles = []
x, y = data1['imag']-data1['zmag'], data1['rmag']-data1['imag']
xerr, yerr = np.sqrt(data1['imagerr']**2+data1['zmagerr']**2), np.sqrt(data1['rmagerr']**2+data1['imagerr']**2)
#print(list(set(data1['date'])))
#print(sorted(list(set(data1['date']))))
marker = 'o'
for i in sorted(list(set(data1['date']))):
    print(i)
    indices = np.where(data1['date'] == i)
    if i == '01282014': indices = np.where( (data1['date'] == i) & (data1['person'] == 'CAT4') )
    if i == '03272014': marker='^'
    plt.errorbar(x[indices], y[indices], xerr=xerr[indices], yerr=yerr[indices], linestyle='None', marker=marker, mec='None', alpha=0.5)
    titles.append(i)
    #plt.xlim(0.2, 1.6)
    #plt.ylim(3,0.5)
    #plt.title(i)
    #plt.show()
plt.legend(titles)
ymin, ymax = plt.ylim()
xmin, xmax = plt.xlim()

t1 = Table.read('iz_ri.txt', format='ascii.commented_header')
#order = np.argsort(t1['i-z'])
Xs = t1['i-z']#[order]
means = t1['r-i_mean']#[order]
sigmas = t1['r-i_sigma']#[order]
plt.plot(Xs, means + 3*sigmas, c='r', lw=1)
plt.plot(Xs, means - 3*sigmas, c='r', lw=1)
plt.plot(Xs, means, 'r:', lw=1)

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xlabel(r'$i-z$')
plt.ylabel(r'$r-i$')
plt.show()
sys.exit()
################

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
