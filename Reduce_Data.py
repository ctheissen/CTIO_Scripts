#!/usr/bin/env python
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
date = '20140321'
biasrange = range(29,46)
gflatrange = range(1,8)
rflatrange = range(8,15)
iflatrange = range(15,22)
zflatrange = range(22,29)
imagerange = range(46,164)

# Things we don't need to change
top = '/mnt/Resources/perseus/CTIO_Data/'

# Try to make a directory for the reduced data
reducedpath = top+date+'/Reduced_Data/'
try: os.mkdir(reducedpath)
except: pass

######################################################################################################### First we make a Master bias
biasData = np.array([fits.getdata(top+date+'/'+date+'.{:0>3}.fits'.format(expId))[:,0:1034] for expId in biasrange])
bias = np.median(biasData, axis=0)
#print(biasData2.shape)
#print(biasData.shape)
#print(bias.shape)

fits.writeto(reducedpath+"MasterBias.fits", bias, clobber=True)

#im = plt.imshow(bias, norm=LogNorm())
#plt.colorbar(im)
#plt.show()
######################################################################################################### Now we do the flats 

# Read in the Master Bias
MasterBias = fits.getdata(reducedpath+'MasterBias.fits')

for band in ['g', 'r', 'i', 'z']: # do all the bands

    flatList = []
    
    # Read the raw flat field exposure from disk
    if band == 'g': ran = gflatrange
    elif band == 'r': ran = rflatrange
    elif band == 'i': ran = iflatrange
    elif band == 'z': ran = zflatrange
            
    for expId in ran:
        
        # Grab the file
        flatRaw = fits.getdata(top+date+'/'+date+'.{:0>3}.fits'.format(expId))[:,0:1034]

        ###### Perform some stats to find the bad pixels. 
        ###### Try clipping the data using the median
        flatRaw2 = stats.sigma_clip(flatRaw, 10)
            
        # Save the masked array to a file so we can read it back later. Only do this on the first flat image
        if expId == 1: flatRaw2.dump(reducedpath+'bad_pixels.pkl')
    
        # The flatfield exposures need to be debiased, just like the data exposure
        dbFlat = flatRaw - MasterBias
            
        # Append our debiased flat field to the list
        flatList.append(dbFlat)
        
    flatList2 = np.array(flatList)

    # Take the median of the flats
    flat = np.median(flatList2, axis=0)

    # Should mask out bad pixels at some point!!!!

    normFlat = flat / np.median(flat.flatten())
    
    fits.writeto(reducedpath+"Master_{0}_flat.fits".format(band), normFlat, clobber=True)
    
    #plt.imshow(flat[:,0:1034], norm=LogNorm())
    #plt.imshow(normFlat[:,0:1034], norm=LogNorm())
    #plt.colorbar()
    #plt.hist(normFlat[:,0:1034].flatten(), bins=np.sqrt(len(normFlat[:,0:1034].flatten())), log=True, histtype='step')
    #plt.show()
    #sys.exit()

######################################################################################################### Now do the rest
# First we need to figure out how many images reside in the directory

AllFiles = glob.glob(top+date+'/*.fits')
length = len(AllFiles)

# Read in the Master Flats
MasterBias = fits.getdata(reducedpath+'MasterBias.fits')
Master_gFlat = fits.getdata(reducedpath+'Master_g_flat.fits')
Master_rFlat = fits.getdata(reducedpath+'Master_r_flat.fits')
Master_iFlat = fits.getdata(reducedpath+'Master_i_flat.fits')
Master_zFlat = fits.getdata(reducedpath+'Master_z_flat.fits')

for i in imagerange: # Now we iterate through the rest of the images and apply corrections

    print('Doing image: %s'%i)
    #if i != 207: continue

    try: hdr = fits.getheader(top+date+'/'+date+'.{:0>3}.fits'.format(i))
    except: continue

    # Get the science data
    scienceRaw = fits.getdata(top+date+'/'+date+'.{:0>3}.fits'.format(i))

    # Read the header to find out what band the image is of, we also want to keep the info
    band = hdr['FILTER2']
        
    # Change the header object
    #print(hdr['OBJECT'])
    hdrobj = hdr['OBJECT']
    hdrobj2 = hdrobj+', REDUCED'
    #print(hdrobj2)
    hdr['OBJECT'] = hdrobj2

    # Now we reduce the data
    if band == 'g': MasterFlat = Master_gFlat
    elif band == 'r': MasterFlat = Master_rFlat
    elif band == 'i': MasterFlat = Master_iFlat
    elif band == 'z': MasterFlat = Master_zFlat
        
    scienceRed = (scienceRaw[:,0:1034] - MasterBias) / MasterFlat
    #print(scienceRaw.shape)
    #print(scienceRed.shape)

    fits.writeto(reducedpath+date+'.{:0>3}.reduced.fits'.format(i), scienceRed, hdr, clobber=True)
        
        
