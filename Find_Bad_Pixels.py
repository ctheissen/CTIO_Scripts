#!/usr/local/bin/env python3
import numpy as np
import sys, os, os.path, time, gc, glob
from astropy.table import Table
from astropy.io import ascii
import photutils as ph
import astropy.io.fits as fits
import astropy.stats as stats
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from matplotlib.colors import LogNorm

# Start timing counter to see how fast it runs
start = time.time()
gc.collect()

date = '20140128'
top = '/mnt/Resources/perseus/CTIO_Data/'
path = top+date+'/'
#reducedpath = top+date+'/Reduced_Data/'

#image = 89

######################################################################################################### Grab the reduced image
Values = []
BadPixels = []
for i in range(1, 29):
    imageData = fits.getdata(path+date+'.f{:0>3}.fits'.format(i))[:,0:1034]
    hdr = fits.getheader(path+date+'.f{:0>3}.fits'.format(i))

    if i == 1:

        ################ Try clipping the data using the median
        data2 = stats.sigma_clip(imageData, 10)
        print(data2.mask)
        
        mask = np.where(imageData < 10000)
        maskedarr = np.zeros(imageData.shape, dtype=bool)
        maskedarr[mask] = 1
        
        plt.figure(101)
        bins = np.logspace(3, 5, 1000)
        plt.hist(imageData.flatten(), bins=bins, histtype='step', log=True, color='r')
        plt.hist(data2.compressed(), bins=bins, histtype='step', log=True, color='b')
        plt.xscale('log')
        #plt.hist(imageData.flatten(), bins=np.sqrt(len(imageData.flatten())), histtype='step', log=True)
        #plt.show()
        #sys.exit()
        
        Values.append(imageData.flatten())
        # Find bad pixels
        imageData2 = imageData.flatten()
        mask2 = np.where(imageData2 < 10000)
        BadPixels.append(mask2[0])
        
    else:
        mask = np.where(imageData < 10000)
        maskedarr2 = np.zeros(imageData.shape, dtype=bool)
        maskedarr2[mask] = 1
        print(np.sum(maskedarr == maskedarr2), len(maskedarr.flatten()), len(maskedarr2.flatten()))
        maskedarr = maskedarr2
        
        Values.append(imageData.flatten())
        imageData2 = imageData.flatten()
        mask2 = np.where(imageData2 < 10000)
        BadPixels.append(mask2[0])
    
    fig = plt.figure(10, figsize=(13,13))
    ax = fig.add_subplot(111)
    ax.imshow(imageData, cmap='Greys', origin='lower', norm=LogNorm())

    fig = plt.figure(11, figsize=(13,13))
    ax = fig.add_subplot(111)
    ax.imshow(data2, cmap='Greys', origin='lower', norm=LogNorm())
    #plt.show()
    #sys.exit()

# Find the common bad pixel values, truly bad pixels should show up in all images
result = set(BadPixels[0])
for s in BadPixels[1:]:
    result.intersection_update(s)

# Put back into a boolean array
Pixels = np.zeros(len(imageData.flatten()), dtype=np.dtype(bool)) # initialize an empty boolean array
Pixels[list(result)] = 1  # set the bad pixel values
Pixels2 = Pixels.reshape(imageData.shape) # reshape array into shape of images

# Dump the mask into a file for future use
Pixels2.dump(top+'bad_pixels.pkl')

# That should be the last step, everything else is pointless
sys.exit()

# Plot all the pixel values
Values = np.array(Values).flatten()
plt.close('all')
bins = np.logspace(3, 5, 1000)
fig = plt.figure(10, figsize=(13,5))
ax = fig.add_subplot(111)
ax.hist(Values, bins=bins, histtype='step', log=True, color='r')
plt.show()
sys.exit()
    
sys.exit()
mean, median, std = stats.sigma_clipped_stats(imageData, sigma=3.0)
print(mean, median, std)
sources = ph.irafstarfind(imageData - median, threshold=5.*std, fwhm=3.0)#, exclude_border=True) 
print(sources)
Xs = sources['xcentroid'][np.where(sources['fwhm']<2)]
Ys = sources['ycentroid'][np.where(sources['fwhm']<2)]

positions = (Xs.data, Ys.data)
apertures = ph.CircularAperture(positions, r=8.)

fig = plt.figure(10, figsize=(13,13))
ax = fig.add_subplot(111)

def find_index_of_nearest_xy(x_array, y_array, x_point, y_point):
    distance = (y_array-y_point)**2 + (x_array-x_point)**2
    print('DISTANCE:', distance, distance.min())
    index = np.where(distance==distance.min())
    return index[0][0]

def onclick(event):
    xcheck, ycheck = event.xdata, event.ydata
    N1, N2 = event.xdata, event.ydata
    mutable['x'] = N1
    mutable['y'] = N2
    plt.close()
    print('onclick points:', xcheck, ycheck)
    
mutable = {}
fig.canvas.mpl_connect('button_press_event', onclick)

norm = ImageNormalize(stretch=SqrtStretch())
ax.imshow(imageData, cmap='Greys', origin='lower', norm=LogNorm())
apertures.plot(color='blue', lw=1.5, alpha=0.5)
plt.show()

# Now let's match up our selection
index1 = find_index_of_nearest_xy(Xs, Ys, mutable['x'], mutable['y'])
"""
fig = plt.figure(10, figsize=(13,13))
ax = fig.add_subplot(111)
apertures = ph.CircularAperture((Xs[index1],Ys[index1]), r=10.)
norm = ImageNormalize(stretch=SqrtStretch())
ax.imshow(imageData, cmap='Greys', origin='lower', norm=LogNorm())
apertures.plot(color='red', lw=1.5, alpha=0.5)
plt.show()
"""

######################################################################################################### Now we do the photomery 

#print(ph.psf.create_prf(imageData, positions, 3))
#Gauss = ph.psf.GaussianPSF(sigma = sources['fwhm'][index1]/2.355s, x_0=Xs[index1], y_0=Ys[index1])

# Figure out how big the aperture should be
for i in range(1,20):
    apertures = ph.CircularAperture((Xs[index1], Ys[index1]), r=i)
    phot_table1 = ph.aperture_photometry(imageData, apertures)
    flux1 = phot_table1['aperture_sum']
    apertures = ph.CircularAperture((Xs[index1], Ys[index1]), r=i+1)
    phot_table2 = ph.aperture_photometry(imageData, apertures)
    flux2 = phot_table2['aperture_sum']
    if 100 - (flux1/flux2 *100) < 1:
        radius = float(i)
        break

print(radius)
"""    
plt.figure(123)
plt.plot(apes, fluxes)
plt.show()
"""

fig = plt.figure(10, figsize=(13,13))
ax = fig.add_subplot(111)
apertures = ph.CircularAperture((Xs[index1],Ys[index1]), r=radius)
norm = ImageNormalize(stretch=SqrtStretch())
ax.imshow(imageData, cmap='Greys', origin='lower', norm=LogNorm())
apertures.plot(color='red', lw=1.5, alpha=0.5)
plt.show()


apertures = ph.CircularAperture((Xs[index1], Ys[index1]), r=radius)
phot_table = ph.aperture_photometry(imageData - median, apertures)
print(phot_table)

# Compute the instrumental magnitude
inM = -2.5*np.log10(phot_table['aperture_sum'] / hdr['EXPTIME'])
print(inM)
