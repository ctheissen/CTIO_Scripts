#!/usr/bin/env python3
import numpy as np
import sys, os, os.path, time, gc, glob, random
from astropy.table import Table
from astropy.io import ascii
import photutils as ph
import astropy.io.fits as fits
import astropy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Start timing counter to see how fast it runs
start = time.time()
gc.collect()

date = '20140128'
top = '/mnt/Resources/perseus/CTIO_Data/'
path = top+date+'/'
reducedpath = top+date+'/Reduced_Data/'

plotT = False

######################################################################################################### Grab the bad pixel mask
masked_array = np.load(reducedpath+'bad_pixels.pkl')
mask = masked_array.mask
######################################################################################################### Grab two bias images

# Going to put this into a function and do it a lot of times

def BiasLoop():

    # Randomly pick two different bias images
    i1 = random.randrange(29,46)
    while True:
        i2 = random.randrange(29,46)
        if i2 != i1: break
    #print('Image Numbers:', i1, i2)

    Data1 = fits.getdata(path+date+'.f{:0>3}.fits'.format(i1))[:,0:1034]
    Data2 = fits.getdata(path+date+'.f{:0>3}.fits'.format(i2))[:,0:1034]
    #hdr = fits.getheader(path+date+'.f{:0>3}.fits'.format(i))

    Data3 = Data2 - Data1
    Data4 = np.ma.array(Data3, mask=mask)
    Data5 = stats.sigma_clip(Data4.compressed(), 3)

    mean, std = np.mean(Data4.compressed()), np.std(Data4.compressed())
    #print('mean','sigma', 'RN(ADU)')
    #print(mean,std, std/np.sqrt(2))
    mean, std = np.mean(Data5.compressed()), np.std(Data5.compressed())
    #print(mean,std, std/np.sqrt(2), '3-sigma clipped')
    """
    # Try to fit a Gaussian
    from scipy.optimize import curve_fit
    from scipy import asarray as ar,exp
    x = ar(range(len(Data5.compressed())))
    y = Data5.compressed()
    def gaus(x,a,x0,sigma): return a*exp(-(x-x0)**2/(2*sigma**2))
    popt,pcov = curve_fit(gaus,x,y,p0=[1,mean,std])
    print(popt)
    """
    if plotT == True:
        plt.figure(10)
        bins = 30# np.sqrt(len(Data4.compressed()))
        plt.hist(Data5.compressed(), bins = bins, histtype='step', log=True)#, range=(-20,20)
        plt.axvline(0, ls=':', c='k')
        plt.show()
    
    return std/np.sqrt(2)

stds = []
for i in range(1000):
    print(i)
    stds.append(BiasLoop())

print(np.mean(stds), np.std(stds)/np.sqrt(len(stds)))
plt.figure(101)
bins = np.sqrt(len(stds))
plt.hist(stds, bins = bins, histtype='step')#, range=(-20,20)
plt.axvline(0, ls=':', c='k')
plt.show()
