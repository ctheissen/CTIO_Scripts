#!/usr/local/bin/env python3
import numpy as np
import sys, os, os.path, time, gc, glob
from astropy.table import Table
import astropy.stats as stats
from astropy.io import ascii
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# WE CAN MAKE THIS SCRIPT BETTER BY PUTTING ALL THE FILES IN AN ARRAY AND THEN LOOKING AT HEADER INFO FOR EACH FILE TO DO EACH STEP

# Things we need to manually set
date       = '20140128'
biasrange  = range(294,305)
gflatrange = range(1,8)
rflatrange = range(8,15)
iflatrange = range(15,22)
zflatrange = range(22,29)
imagerange = range(29,294)

# Things we don't typically need to change
top = '/mnt/Resources/perseus/CTIO_Data/'

# Try to make a directory for the reduced data
reducedpath = top+date+'/Reduced_Data/'
try: os.mkdir(reducedpath)
except: pass

######################################################################################################### First we make a Master bias
# Try to grab the file
Passed = 0 # Set identifier to see if we found the file
try:
    biasData = np.array([fits.getdata(top+date+'/'+date+'.{:0>3}.fits'.format(expId))[:,0:1034] for expId in biasrange])
    Passed = 1
except:
    print('Trying a different file name.')
    Passed = 0
    pass
if Passed == 0: # We did not find the file the first time.
    try: biasData = np.array([fits.getdata(top+date+'/'+date+'.f{:0>3}.fits'.format(expId))[:,0:1034] for expId in biasrange])
    except: raise OSError('We do not know that filename: %s'%(top+date+'/'+date+'.f{:0>3}.fits'.format(expId)))

# Take the median across the cube
bias = np.median(biasData, axis=0)

# Write the master
fits.writeto(reducedpath+"MasterBias.fits", bias, clobber=True)
######################################################################################################### Now we do the flats 

# Read in the Master Bias
MasterBias = fits.getdata(reducedpath+'MasterBias.fits')

for band in ['g', 'r', 'i', 'z']: # do all the bands

    # Make an empty list to put the arrays in
    flatList = []
    
    # Read the raw flat field exposure from disk
    if band == 'g': ran = gflatrange
    elif band == 'r': ran = rflatrange
    elif band == 'i': ran = iflatrange
    elif band == 'z': ran = zflatrange
            
    for expId in ran:
        
        # Try to grab the file
        Passed = 0 # Set identifier to see if we found the file
        try:
            flatRaw = fits.getdata(top+date+'/'+date+'.{:0>3}.fits'.format(expId))[:,0:1034] # we are omitting bad pixels
            Passed = 1
        except:
            print('Trying a different file name.')
            Passed = 0
            pass
        if Passed == 0: # We did not find the file the first time.
            try: flatRaw = fits.getdata(top+date+'/'+date+'.f{:0>3}.fits'.format(expId))[:,0:1034] # we are omitting bad pixels
            except: raise OSError('We do not know that filename: %s'%(top+date+'/'+date+'.f{:0>3}.fits'.format(expId)))
            
        # The flatfield exposures need to be debiased, just like the data exposure
        dbFlat = flatRaw - MasterBias
            
        # Append our debiased flat field to the list
        flatList.append(dbFlat)

    # Convert to a Numpy array
    flatList2 = np.array(flatList)

    # Take the median of the flats
    flat = np.median(flatList2, axis=0)

    # Normalize
    normFlat = flat / np.median(flat.flatten())

    # Write the master flat
    fits.writeto(reducedpath+"Master_{0}_flat.fits".format(band), normFlat, clobber=True)

######################################################################################################### Now do the rest
# First we need to figure out how many images reside in the directory

AllFiles = glob.glob(top+date+'/*.fits')
length   = len(AllFiles)

# Read in the Master Flats
MasterBias   = fits.getdata(reducedpath+'MasterBias.fits')
Master_gFlat = fits.getdata(reducedpath+'Master_g_flat.fits')
Master_rFlat = fits.getdata(reducedpath+'Master_r_flat.fits')
Master_iFlat = fits.getdata(reducedpath+'Master_i_flat.fits')
Master_zFlat = fits.getdata(reducedpath+'Master_z_flat.fits')

for i in imagerange: # Now we iterate through the rest of the images and apply corrections

    print('Doing image: %s'%i)
    #if i != 207: continue

    # Try to grab the file
    Passed = 0
    try:
        hdr = fits.getheader(top+date+'/'+date+'.{:0>3}.fits'.format(i))
        Passed = 1
    except:
        print('Trying a different file name.')
        Passed = 0
        pass
    if Passed == 0:
        try: hdr = fits.getheader(top+date+'/'+date+'.f{:0>3}.fits'.format(i))
        except: raise OSError('We do not know that filename: %s'%(top+date+'/'+date+'.f{:0>3}.fits'.format(i)))
        
    # Try to grab the file
    # Get the science data
    Passed = 0
    try:
        scienceRaw = fits.getdata(top+date+'/'+date+'.{:0>3}.fits'.format(i))
        Passed = 1
    except: 
        print('Trying a different file name.')
        Passed = 0
        pass
    if Passed == 0:
        try: scienceRaw = fits.getdata(top+date+'/'+date+'.f{:0>3}.fits'.format(i))
        except: raise OSError('We do not know that filename: %s'%(top+date+'/'+date+'.f{:0>3}.fits'.format(i)))
        
    # Read the header to find out what band the image is of, we also want to keep the info
    band = hdr['FILTER2']
        
    # Change the header object
    hdrobj = hdr['OBJECT']
    hdrobj2 = hdrobj+', REDUCED'
    hdr['OBJECT'] = hdrobj2

    # Now we reduce the data
    if band == 'g': MasterFlat = Master_gFlat
    elif band == 'r': MasterFlat = Master_rFlat
    elif band == 'i': MasterFlat = Master_iFlat
    elif band == 'z': MasterFlat = Master_zFlat

    # Reduce the raw science data
    scienceRed = (scienceRaw[:,0:1034] - MasterBias) / MasterFlat

    # Write the reduced data to file
    fits.writeto(reducedpath+date+'.{:0>3}.reduced.fits'.format(i), scienceRed, hdr, clobber=True)
        
        
