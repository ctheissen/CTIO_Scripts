#!/usr/bin/env python
import numpy as np
import sys, os, os.path, time, gc, glob
from astropy.table import Table
from astropy.io import ascii
import photutils as ph
import astropy.io.fits as fits
import astropy.stats as stats
import matplotlib as mat
mat.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.wcs import WCS
from pylab import get_current_fig_manager

# Things we need to input
date        = '20140321'  # Night of the observations
image       = 163          # Which object are we doing?
airmass     = 1.301       # What is the airmass?

# Things we don't need to input
top         = '/mnt/Resources/perseus/CTIO_Data/'
reducedpath = top+date+'/Reduced_Data/'

# Tweak params for source extraction if we need to
thres  = 10.   # 5. usually
fwhm   = 3.0   # 3.0 usually

# Other params
Gain = 2.7 # e/ADU
RN   = 1.6 # ADU, or 4.0 electrons (this is per pixel). We empirically measured this to be 1.55.

######################################################################################################### Some functions for use 
def find_index_of_nearest_xy(x_array, y_array, x_point, y_point):
    distance = (y_array-y_point)**2 + (x_array-x_point)**2
    #print('DISTANCE:', distance, distance.min())
    index = np.where(distance==distance.min())
    return index[0][0]

def onclick(event):
    xcheck, ycheck = event.xdata, event.ydata
    N1, N2 = event.xdata, event.ydata
    mutable['x'] = N1
    mutable['y'] = N2
    plt.close("all")
    print('onclick points:', xcheck, ycheck)
    
def onclickclose(event):
    plt.close("all")
    
######################################################################################################### Function to grab the reduced image
def GetImage(image, mask):
    imageData = fits.getdata(reducedpath+date+'.{:0>3}.reduced.fits'.format(image))
    hdr = fits.getheader(reducedpath+date+'.{:0>3}.reduced.fits'.format(image))
    w = WCS(reducedpath+date+'.{:0>3}.reduced.fits'.format(image))

    # Parse the header to get the object name
    ObjName = hdr['OBJECT'].split(',')[0].split(' ')[0]
    band = hdr['FILTER2']

    # Create the masked Image Data
    imageDataM = np.ma.array(imageData, mask=mask)

    # Computed the background levels
    mean, median, std = stats.sigma_clipped_stats(imageData, mask=mask, sigma=3.0)
    print('mean', 'median', 'std', 'BACKGROUND')
    print(mean, median, std)

    # Remove the background
    imageDataRed = imageDataM - median
    return imageDataRed, median, ObjName, band, std, hdr, w
######################################################################################################### Function to do photometry
def doPhot(image, median, Xs, Ys, hdr):
    # Maybe implement PSF photometry sometime in the future
    #print(ph.psf.create_prf(imageData, positions, 3))
    #Gauss = ph.psf.GaussianPSF(sigma = sources['fwhm'][index1]/2.355s, x_0=Xs[index1], y_0=Ys[index1])

    # Fix the masked values so they are just 0. This is not working in photutils,
    image.filled(0)

    # Figure out how big the aperture should be (when the change is < 1% we can stop)
    # WE CAN CLEAN THIS UP BY DOING MULTIPLE APERTURES AT ONCE AND COMPARING AT THE END, AS AN ENSEMBLE
    radius = None
    for i in range(1,20):
        apertures = ph.CircularAperture((Xs, Ys), r=i)
        phot_table1 = ph.aperture_photometry(image, apertures)#, mask=mask)
        flux1 = phot_table1['aperture_sum']
        apertures = ph.CircularAperture((Xs, Ys), r=i+1)
        phot_table2 = ph.aperture_photometry(image, apertures)#, mask=mask)
        flux2 = phot_table2['aperture_sum']
        if 100 - (flux1/flux2 *100) < 1:
            radius = float(i)
            break
    if radius == None: radius=10 # This is the optimal aperture size historically. This is just a fix for contaminted areas.

    fig = plt.figure(10, figsize=(12,12))
    ax = fig.add_subplot(111)
    apertures = ph.CircularAperture((Xs,Ys), r=radius)
    fig.canvas.mpl_connect('button_press_event', onclickclose)
    ax.imshow(image, cmap='Greys', norm=LogNorm())
    apertures.plot(color='red', lw=1.5, alpha=0.5)
    plt.show()

    apertures = ph.CircularAperture((Xs, Ys), r=radius)
    phot_table = ph.aperture_photometry(image, apertures)#, mask=mask)
    #print(phot_table)

    # Compute the instrumental magnitude and associated error
    inM = -2.5*np.log10(phot_table['aperture_sum'].data / hdr['EXPTIME'])
    sigma_mag = 1.0857 * np.sqrt( phot_table['aperture_sum'].data*Gain + 2*np.pi*radius**2 * (median1*Gain + (RN*Gain)**2) ) / (phot_table['aperture_sum'].data*Gain)
    return inM, sigma_mag
######################################################################################################### Grab the bad pixel mask
#masked_array = np.load(top+'bad_pixels.pkl')
#mask = masked_array.mask
mask = np.load(top+'bad_pixels.pkl')
######################################################################################################### Grab the corrections
Corr = ascii.read(reducedpath+'Corrections.txt')
######################################################################################################### Grab the list of objects
AllFiles = glob.glob(reducedpath+date+'*.fits') # Get a list of all the relevant file names
FileNums = []
Objs = []
Bands = []
for File in AllFiles:
    
    # Get the file number
    FileNum = int(File.split('/')[-1].split('.')[1].strip('f'))
    FileNums.append(FileNum)
    
    # Get the object name
    hdr = fits.getheader(File)
    Obj = hdr['OBJECT'].split(',')[0].split(' ')[0]
    Objs.append(Obj)

    # Get the band
    Bands.append(hdr['FILTER2'])

# Convert to numpy arrays for use later
Objs = np.array(Objs)
FileNums = np.array(FileNums)
Bands = np.array(Bands)
######################################################################################################### Grab the reduced image

# First let's get the images and the necessary data associated with them
image1, median1, ObjName1, band1, std1, hdr1, w1 = GetImage(image, mask)
print(hdr1['OBJECT'])
Object1 = hdr1['OBJECT'].split(',')[0].split(' ')[0] 

##################################### Grab the first and last image to compare brightness
# Get the images with our object
Relevant = np.where(Objs == Object1)

# Get the images with our object
firstImage = FileNums[Relevant][0]

# Figure out what the reddest band is
BandCheck = 'g'
Bandz = np.where(Bands[Relevant] == 'z')
if len(Bandz) == 0:
    Bandi = np.where(Bands[Relevant] == 'i')
    BandCheck = 'i'
    if len(Bandi) == 0:
        Bandr = np.where(Bands[Relevant] == 'r')
        BandCheck = 'r'
else: BandCheck = 'z'

if BandCheck == 'z': lastImage = FileNums[Relevant][np.where(Bands[Relevant] == 'z')][-1]
elif BandCheck == 'i': lastImage = FileNums[Relevant][np.where(Bands[Relevant] == 'i')][-1]
elif BandCheck == 'r': lastImage = FileNums[Relevant][np.where(Bands[Relevant] == 'r')][-1]

imageShow1, medianShow1, ObjNameShow1, bandShow1, stdShow1, hdrShow1, wShow1 = GetImage(firstImage, mask)
imageShow2, medianShow2, ObjNameShow2, bandShow2, stdShow2, hdrShow2, wShow2 = GetImage(lastImage, mask)

# Plot the bluest and reddest image for visual comparison
figLook = plt.figure(187, figsize=(10,5))
#figLook.canvas.mpl_connect('button_press_event', onclickclose)
ax1 = figLook.add_subplot(121)
ax2 = figLook.add_subplot(122)
ax1.imshow(imageShow1, cmap='Greys', norm=LogNorm())
ax1.set_title('%s-band'%hdrShow1['FILTER2'])
ax2.imshow(imageShow2, cmap='Greys', norm=LogNorm())
ax2.set_title('%s-band'%hdrShow2['FILTER2'])

# Move the figure on the screen
thismanager = get_current_fig_manager()
thismanager.window.wm_geometry("+1100+0")
#figLook.canvas.manager.window.SetPosition((500, 0))
#####################################

# Choose how to do source extraction
threshold = thres*std1
#sources = ph.irafstarfind(imageDataRed, threshold=threshold, fwhm=fwhm)#, exclude_border=True) 
sources = ph.daofind(image1, threshold=threshold, fwhm=fwhm)#, exclude_border=True) 
#print(sources)

# Put a mask here that we might not need anymore
Xs = sources['xcentroid']#[np.where(sources['fwhm']<2)]
Ys = sources['ycentroid']#[np.where(sources['fwhm']<2)]

positions = (Xs.data, Ys.data)
apertures = ph.CircularAperture(positions, r=10.)

fig = plt.figure(10, figsize=(12,12))
ax = fig.add_subplot(111)

mutable = {}
fig.canvas.mpl_connect('button_press_event', onclick)

ax.imshow(image1, cmap='Greys', norm=LogNorm())
apertures.plot(color='blue', lw=1.5, alpha=0.5)
plt.show()

# Now let's match up our selection
index1 = find_index_of_nearest_xy(Xs, Ys, mutable['x'], mutable['y'])

######################################################################################################### Now we do the photomery 

Mag, unMag = doPhot(image1, median1, Xs[index1], Ys[index1], hdr1)

# Compute the final magnitude and uncertainty
if hdr1['FILTER2'] == 'g': 
    TrueMag   = Corr['deltaM'][0]+Mag
    unTrueMag = np.sqrt(unMag**2 + Corr['deltaMun'][0]**2 + airmass**2*Corr['Kun'][0]**2)
elif hdr1['FILTER2'] == 'r': 
    TrueMag   = Corr['deltaM'][1]+Mag
    unTrueMag = np.sqrt(unMag**2 + Corr['deltaMun'][1]**2 + airmass**2*Corr['Kun'][1]**2)
elif hdr1['FILTER2'] == 'i': 
    TrueMag   = Corr['deltaM'][2]+Mag
    unTrueMag = np.sqrt(unMag**2 + Corr['deltaMun'][2]**2 + airmass**2*Corr['Kun'][2]**2)
elif hdr1['FILTER2'] == 'z': 
    TrueMag   = Corr['deltaM'][3]+Mag
    unTrueMag = np.sqrt(unMag**2 + Corr['deltaMun'][3]**2 + airmass**2*Corr['Kun'][3]**2)

print(Object1, hdr1['FILTER2'], TrueMag[0], unTrueMag[0])

# Check if file exists and write some stuff
if os.path.isfile(reducedpath+'Sources.txt'):
    f = open(reducedpath+'Sources.txt', 'a')
    f.write('%s,%s,%s,%s,\n'%(Object1, hdr1['FILTER2'], TrueMag[0], unTrueMag[0]))
    f.close()
else:
    f = open(reducedpath+'Sources.txt', 'a')
    f.write('#Object,Band,Mag,MagUn,Notes\n')
    f.write('%s,%s,%s,%s,\n'%(Object1, hdr1['FILTER2'], TrueMag[0], unTrueMag[0]))
    f.close()
