#!/usr/bin/env python
import numpy as np
import scipy as sp
import sys, os, os.path, time, gc, glob
from astropy.table import Table, hstack, vstack
from astropy.io import ascii
import photutils as ph
import astropy.io.fits as fits
import astropy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pylab as py
from astropy.wcs import WCS
import matplotlib.figure as mfigure
from matplotlib.backends.backend_agg import FigureCanvasAgg #canvas

# Things we need to manually input
date      = '20140321'               # Night we are doing
start     = [60, 86, 99]           # Starting (g-band) image numbers
plus      = 3                        # 0 for g-band, 1 for r-band, 2 for i-band, 3 for z-band
airmasses = [1.119, 1.292, 1.652]    # Need to manually input this since it's not in the headers

# Things we don't need to change
top         = '/mnt/Resources/perseus/CTIO_Data/' # path to CTIO data
reducedpath = top+date+'/Reduced_Data/'           # path to reduced CTIO data
image0    = start[0]+plus                         # starting image for first set of calibrations
image00   = start[1]+plus                         # starting image for second set of calibrations
image000  = start[2]+plus                         # starting image for third set of calibrations

# Tweak params for source extraction if we need to
thres = 5.  # 5. usually
fwhm = 3.0  # 3.0 usually

# Other params that are used
Gain = 2.7 # e/ADU
RN = 1.6 # ADU, or 4.0 electrons (this is per pixel). We empirically measured this to be 1.55.

######################################################################################################### Functions for use
def overlap(a,b):
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]

def find_index_of_nearest_xy(x_array, y_array, x_point, y_point):
    distance = (y_array-y_point)**2 + (x_array-x_point)**2
    #print('DISTANCE:', distance, distance.min())
    index = np.where(distance==distance.min())
    return index[0][0]

def find_index_of_nearest_xy2(x_array, y_array, x_point, y_point):
    distance = (y_array-y_point)**2 + (x_array-x_point)**2
    #print('DISTANCE:', np.sqrt(distance), np.sqrt(distance.min()))
    index = np.where(distance==distance.min())
    if np.sqrt(distance.min()) < 4: return index[0][0]
    else: return

def find_index_of_nearest_radec(ra_array, dec_array, ra_point, dec_point):
    distance = (dec_array-dec_point)**2 + (ra_array-ra_point)**2
    #print('DISTANCE:', np.sqrt(distance), np.sqrt(distance.min()))
    #print('MIN DISTANCE:', np.sqrt(distance.min()))
    index = np.where(distance==distance.min())
    if np.sqrt(distance.min()) < 3e-4: return index[0][0]
    else: return 

def clicky(event):
    if event.key in ('e', 'E','q','Q'): # Kill it when a quit key is pressed
        plt.close()
        return
    XsM.append(event.xdata)
    YsM.append(event.ydata)
    ax.scatter(event.xdata, event.ydata, marker='x', c='r', s=60)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    plt.draw()
    print('onclick points:', event.xdata, event.ydata)

def onclickclose(event):
    plt.close()

######################################################################################################### Function to grab the reduced image
def GetImage(image, mask):
    imageData = fits.getdata(reducedpath+date+'.{:0>3}.reduced.new'.format(image))
    hdr = fits.getheader(reducedpath+date+'.{:0>3}.reduced.new'.format(image))
    w = WCS(reducedpath+date+'.{:0>3}.reduced.new'.format(image))

    # Parse the header to get the object name
    ObjName = hdr['OBJECT'].split(',')[0].split(' ')[0]
    band = hdr['FILTER2']

    # Create the masked Image Data
    imageDataM = np.ma.array(imageData, mask=mask)

    # Computed the background levels
    mean, median, std = stats.sigma_clipped_stats(imageData, mask=mask, sigma=3.0)
    print('mean', 'median', 'std', 'BACKGROUND')
    print(mean, median, std)

    ## Remove the background
    imageDataRed = imageDataM - median
    #imageDataRed = imageDataM # or do not

    #return imageDataRed, median, ObjName, band, std, hdr, w
    return imageDataRed, median, ObjName, band, std, hdr, w

######################################################################################################### Function to match coords across images
def match2(RAs1, DECs1, RAs2, DECs2):
    indices0 = []
    for ra1, dec1 in zip(RAs1, DECs1):
        index1 = find_index_of_nearest_radec(RAs2, DECs2, ra1, dec1)
        if index1 != None: indices0.append(index1)
    indices1 = np.array(indices0)
    return indices1

def match2sdss(RAs1, DECs1, RAsSDSS, DECsSDSS):
    indices0 = []
    for ra1, dec1 in zip(RAs1, DECs1):
        index1 = find_index_of_nearest_radec(RAsSDSS, DECsSDSS, ra1, dec1)
        if index1 != None: indices0.append(index1)
    indices1 = np.array(indices0)
    return indices1

def match3(RAs1, DECs1, RAs2, DECs2, RAs3, DECs3):
    indices0 = []
    for ra1, dec1 in zip(RAs2, DECs2):
        index1 = find_index_of_nearest_radec(RAs1, DECs1, ra1, dec1)
        if index1 != None: indices0.append(index1)
    indices00 = []
    for ra1, dec1 in zip(RAs3, DECs3):
        index1 = find_index_of_nearest_radec(RAs1, DECs1, ra1, dec1)
        if index1 != None: indices00.append(index1)
    overlap1 = np.nonzero(np.in1d(indices0, indices00))[0]
    indices1 = np.array(indices0)[overlap1]
    return indices1

######################################################################################################### Function to do photometry
def doPhot(image, median, Xs, Ys, hdr):
    # Maybe implement PSF photometry sometime in the future
    #print(ph.psf.create_prf(imageData, positions, 3))
    #Gauss = ph.psf.GaussianPSF(sigma = sources['fwhm'][index1]/2.355s, x_0=Xs[index1], y_0=Ys[index1])

    # Fix the masked values so they are just 0. This is not working in photutils,
    image.filled(0)

    # Figure out how big the aperture should be (when the change is < 1% we can stop)
    radii = []
    for j in range(len(Xs)): # Loop over every extracted source
        for i in range(1,20):
            apertures = ph.CircularAperture((Xs[j], Ys[j]), r=i)
            phot_table1 = ph.aperture_photometry(image, apertures)#, mask=mask)
            flux1 = phot_table1['aperture_sum']
            apertures = ph.CircularAperture((Xs[j], Ys[j]), r=i+1)
            phot_table2 = ph.aperture_photometry(image, apertures)#, mask=mask)
            flux2 = phot_table2['aperture_sum']
            if 100 - (flux1/flux2 *100) < 1:
                radii.append(float(i))
                break
            
    # Now we grab the most common aperture size
    try: radius = sp.stats.mode(radii)[0][0]
    except: radius = 10 # This is the typical optimal value
    print('Radius size:', radius)

    # Show the objects and the aperture size
    fig = plt.figure(10, figsize=(12,12))
    ax = fig.add_subplot(111)
    apertures = ph.CircularAperture((Xs,Ys), r=radius)
    fig.canvas.mpl_connect('button_press_event', onclickclose)
    ax.imshow(image, cmap='Greys', norm=LogNorm())
    apertures.plot(color='red', lw=1.5, alpha=0.5)
    plt.show()

    # Now let's do the photometry
    apertures = ph.CircularAperture((Xs, Ys), r=radius) # Create circular apertures for each source
 
    phot_table = ph.aperture_photometry(image, apertures)#, mask=mask)
    #print(phot_table)

    # Compute the instrumental magnitude and associated error
    inM = -2.5*np.log10(phot_table['aperture_sum'].data / hdr['EXPTIME'])
    #print(inM)

    sigma_mag = 1.0857 * np.sqrt( phot_table['aperture_sum'].data*Gain + 2*np.pi*radius**2 * (median*Gain + (RN*Gain)**2) ) / (phot_table['aperture_sum'].data*Gain)
    #print(sigma_mag)
    #print(ObjName2, band2, hdr2['EXPTIME'], inM, sigma_mag)
    return inM, sigma_mag

######################################################################################################### Grab the bad pixel mask
#masked_array = np.load(top+'bad_pixels.pkl')
#mask = masked_array.mask
mask = np.load(top+'bad_pixels.pkl')
######################################################################################################### Now let's match to the SDSS Calib Field
Field0 = ascii.read(top+'DR10_SDSS_CALIB_FIELD1.csv')
#Field1 = Field0[np.where(Field0['rmag'] < 18)]
Field00 = ascii.read(top+'DR10_SDSS_CALIB_FIELD2.csv')
#Field3 = Field2[np.where(Field2['rmag'] < 18)]
######################################################################################################### Start the show

# First let's get the images and the necessary data associated with them
image1, median1, ObjName1, band1, std1, hdr1, w1 = GetImage(image0, mask)
image2, median2, ObjName2, band2, std2, hdr2, w2 = GetImage(image00, mask)
image3, median3, ObjName3, band3, std3, hdr3, w3 = GetImage(image000, mask)

# Grab sources in the first image
threshold = thres*std1
sources = ph.daofind(image1, threshold=threshold, fwhm=fwhm) 
# Filter out sources close to the edges
sX, sY = sources['xcentroid'][np.where( (sources['xcentroid']>=20) & (sources['xcentroid']<=1000) & (sources['ycentroid']>=20) & (sources['ycentroid']<=1000) )], sources['ycentroid'][np.where( (sources['xcentroid']>=20) & (sources['xcentroid']<=1000) & (sources['ycentroid']>=20) & (sources['ycentroid']<=1000) )]
Xs1 = sX
Ys1 = sY

# Grab sources in the second image
threshold = thres*std2
sources = ph.daofind(image2, threshold=threshold, fwhm=fwhm)
# Filter out sources close to the edges
sX, sY = sources['xcentroid'][np.where( (sources['xcentroid']>=20) & (sources['xcentroid']<=1000) & (sources['ycentroid']>=20) & (sources['ycentroid']<=1000) )], sources['ycentroid'][np.where( (sources['xcentroid']>=20) & (sources['xcentroid']<=1000) & (sources['ycentroid']>=20) & (sources['ycentroid']<=1000) )]
Xs2 = sX
Ys2 = sY

# Grab sources in the third image
threshold = thres*std3
sources = ph.daofind(image3, threshold=threshold, fwhm=fwhm)
# Filter out sources close to the edges
sX, sY = sources['xcentroid'][np.where( (sources['xcentroid']>=20) & (sources['xcentroid']<=1000) & (sources['ycentroid']>=20) & (sources['ycentroid']<=1000) )], sources['ycentroid'][np.where( (sources['xcentroid']>=20) & (sources['xcentroid']<=1000) & (sources['ycentroid']>=20) & (sources['ycentroid']<=1000) )]
Xs3 = sX
Ys3 = sY

#################################### Now let's find sources in common between all the images
# First convert everything to RADEC for comparison
RAs1, DECs1 = w1.wcs_pix2world(Xs1,Ys1, 0)
RAs2, DECs2 = w2.wcs_pix2world(Xs2,Ys2, 0)
RAs3, DECs3 = w3.wcs_pix2world(Xs3,Ys3, 0)

indices1 = match3(RAs1, DECs1, RAs2, DECs2, RAs3, DECs3)
indices2 = match3(RAs2, DECs2, RAs1, DECs1, RAs3, DECs3)
indices3 = match3(RAs3, DECs3, RAs2, DECs2, RAs1, DECs1)
####################################

# Make the apertures
positions1 = (Xs1.data[indices1], Ys1.data[indices1])
apertures1 = ph.CircularAperture(positions1, r=10.)
positions2 = (Xs2.data[indices2], Ys2.data[indices2])
apertures2 = ph.CircularAperture(positions2, r=10.)
positions3 = (Xs3.data[indices3], Ys3.data[indices3])
apertures3 = ph.CircularAperture(positions3, r=10.)

# Put up a plot of each band just for reference
fig = plt.figure(102, figsize=(20,5))
fig.canvas.mpl_connect('button_press_event', onclickclose)
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
ax1.imshow(image1, aspect='equal', cmap='Greys', norm=LogNorm())
ax2.imshow(image2, aspect='equal', cmap='Greys', norm=LogNorm())
ax3.imshow(image3, aspect='equal', cmap='Greys', norm=LogNorm())
apertures1.plot(ax=ax1, color='red', lw=1.5, alpha=0.5)
apertures2.plot(ax=ax2, color='red', lw=1.5, alpha=0.5)
apertures3.plot(ax=ax3, color='red', lw=1.5, alpha=0.5)

################################### Convert SDSS coords to pixel coords and plot on top so we know which stars have overlap (using image1)
###### First, figure out which field to use
# Is there ovelap between the field and (one of) our objects?
RAcheck, DECcheck = w1.wcs_pix2world(Xs1.data[indices1][0], Ys1.data[indices1][0], 0)
C = [Field0['ra'][0]-1, Field0['ra'][0]+1]
M = [RAcheck-.1, RAcheck+.1]
overlap1 = overlap(C, M)
C = [Field0['dec'][0]-1, Field0['dec'][0]+1]
M = [DECcheck-.1, DECcheck+.1]
overlap2 = overlap(C, M)

# Test which field we are looking at
if overlap1 and overlap2: Field000 = Field0
else: Field000 = Field00

Field1 = Field000[np.where( (Field000['gmagerr']< .1) & (Field000['rmagerr']< .1) & (Field000['imagerr']< .1) & (Field000['zmagerr']< .1) )] # Only get good stars
###################################

xSDSS2, ySDSS2 = w1.wcs_world2pix(Field1['ra'], Field1['dec'], 0)
#print(len(Field1), len(xSDSS2))
xSDSS = xSDSS2[np.where( (xSDSS2 > 0) & (xSDSS2 < 1034) & (ySDSS2 > 0) & (ySDSS2 < 1034) )]
ySDSS = ySDSS2[np.where( (xSDSS2 > 0) & (xSDSS2 < 1034) & (ySDSS2 > 0) & (ySDSS2 < 1034) )]

#################################### Now let's find sources in common between all the images and SDSS
# Now let's match up our selections to extracted sources
indicesSDSS = []
for x, y in zip(Xs1.data[indices1], Ys1.data[indices1]):
    if find_index_of_nearest_xy2(xSDSS, ySDSS, x, y) == None: continue
    else: indicesSDSS.append(find_index_of_nearest_xy2(xSDSS, ySDSS, x, y))

print(len(xSDSS), len(ySDSS), len(indicesSDSS))
#print(xSDSS, indicesSDSS)
####################################
XsM = []
YsM = []

fig = plt.figure(10, figsize=(12,12))
ax = fig.add_subplot(111)
ax.imshow(image1, aspect='equal', cmap='Greys', norm=LogNorm())
apertures1.plot(ax=ax, color='red', lw=1.5, alpha=0.5)
# Plot SDSS
positionsSDSS = (xSDSS, ySDSS)
aperturesSDSS = ph.CircularAperture(positionsSDSS, r=10.)
aperturesSDSS.plot(ax=ax, color='blue', lw=1.5, alpha=0.5)
#cid0 = plt.gcf().canvas.mpl_connect('button_press_event', clicky)
#cid1 = plt.gcf().canvas.mpl_connect('key_press_event', clicky)
fig.canvas.mpl_connect('button_press_event', onclickclose)
xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()
ax.set_title('Press q, Q, e, or E to exit after clicking stars')
ax.scatter(xSDSS[indicesSDSS], ySDSS[indicesSDSS], marker='x', c='r', s=60)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
plt.show()

#print(XsM, YsM)
XsM, YsM = xSDSS[indicesSDSS], ySDSS[indicesSDSS]

# Now let's match up our selections to extracted sources
indices = []
for x,y in zip(XsM, YsM):
    indices.append(find_index_of_nearest_xy(Xs1, Ys1, x, y))

#print(indices)
#print(Xs1[indices], Ys1[indices])

#################### Now we need to match it to the other 2 images
Xs11, Ys11 = Xs1[indices], Ys1[indices]
RAs11, DECs11 = w1.wcs_pix2world(Xs11,Ys11, 0)
indices22 = match2(RAs11, DECs11, RAs2, DECs2)
indices33 = match2(RAs11, DECs11, RAs3, DECs3)
Xs22, Ys22 = Xs2[indices22], Ys2[indices22]
Xs33, Ys33 = Xs3[indices33], Ys3[indices33]

######################################################################################################### Now we do the photomery 
Mags1, unMags1 = doPhot(image1, median1, Xs11, Ys11, hdr1)
Mags2, unMags2 = doPhot(image2, median2, Xs22, Ys22, hdr2)
Mags3, unMags3 = doPhot(image3, median3, Xs33, Ys33, hdr3)

#print(Mags1, unMags1)
#print(Mags2, unMags2)
#print(Mags3, unMags3)

# Plot to see the extinction curve
Ks = []
fig = plt.figure(555)
fig.canvas.mpl_connect('button_press_event', onclickclose)
for i in range(len(Mags1)):
    x = airmasses
    y = [Mags1[i], Mags2[i], Mags3[i]]
    yerr = [unMags1[i], unMags2[i], unMags3[i]]
    z = np.polyfit(x, y, 1, w = 1/np.array(yerr)**2)
    print('COEFFS:', z)
    p = np.poly1d(z)
    xp = np.linspace(0, 2)
    plt.errorbar(x, y, yerr = yerr, fmt='o')
    plt.plot(xp, p(xp))
    Ks.append(z[0])
#plt.show()

# Remove the NaNs and make some sigma clipped data
#print('Ks', Ks, len(Ks))
Ks = np.array(Ks)
Ks2 = Ks[~np.isnan(Ks)]
filtered_data0 = stats.sigma_clip(Ks2, 2, None)
filtered_data = filtered_data0[~filtered_data0.mask]

plt.close()
fig = plt.figure(555)
fig.canvas.mpl_connect('button_press_event', onclickclose)
plt.hist(Ks2, histtype='step', color='r')
plt.hist(filtered_data, histtype='step', color='b')
print(np.mean(Ks2), np.median(Ks2), np.std(Ks2), np.std(Ks2)/np.sqrt(len(Ks2)))
print(np.mean(filtered_data), np.median(filtered_data), np.std(filtered_data), np.std(filtered_data)/np.sqrt(len(filtered_data)))
plt.legend(['Unfiltered','Filtered'])
plt.xlabel('Extinction Coefficient')
plt.ylabel('Number')
plt.show()

K = np.mean(filtered_data)
unK = np.std(filtered_data)/np.sqrt(len(filtered_data))
print('Number of stars used for extinction calibration: %s'%len(filtered_data))

# Now we take our K value and extrapolate to zero atmosphere, then compare to SDSS magnitude
Mags0_1 = Mags1[~np.isnan(Ks)][~filtered_data0.mask]  - K*airmasses[0]
unMags0_1 = np.sqrt(unMags1[~np.isnan(Ks)][~filtered_data0.mask]**2 + airmasses[0]**2*unK**2)
#print(Mags0_1, unMags0_1)

##################################################################### SDSS Calibration

# Just take the SDSS stars we care about
Field11 = Field1[np.where( (xSDSS2 > 0) & (xSDSS2 < 1034) & (ySDSS2 > 0) & (ySDSS2 < 1034) )][indicesSDSS][~np.isnan(Ks)][~filtered_data0.mask]
#print(Field11['rmag'])
if band1 == 'g': TrueMags, unTrueMags = Field11['gmag'], Field11['gmagerr']
elif band1 == 'r': TrueMags, unTrueMags = Field11['rmag'], Field11['rmagerr']
elif band1 == 'i': TrueMags, unTrueMags = Field11['imag'], Field11['imagerr']
elif band1 == 'z': TrueMags, unTrueMags = Field11['zmag'], Field11['zmagerr']
deltaMag = TrueMags - Mags0_1

# Fit a straight line to the points
z = np.polyfit(TrueMags, Mags0_1, 1)
print('COEFFS:', z)
p = np.poly1d(z)
xp = np.linspace(TrueMags.min(), TrueMags.max())

plt.close()
fig = plt.figure(555)
fig.canvas.mpl_connect('button_press_event', onclickclose)
plt.errorbar(TrueMags, Mags0_1, xerr = unTrueMags, yerr=unMags0_1, ls='None')
plt.plot(xp, p(xp), 'r--')
plt.xlabel('SDSS Mag')
plt.ylabel('Instrumental Mag (Outside Atmosphere)')
plt.show()


# Filter again? (Not sure I should do this, feels like fudging the numbers)
filtered_data00 = stats.sigma_clip(deltaMag, 2, None)
filtered_data2 = filtered_data00[~filtered_data00.mask]
print(deltaMag)
print(np.mean(deltaMag), np.std(deltaMag)/np.sqrt(len(deltaMag)))
print(np.mean(filtered_data2), np.std(filtered_data2)/np.sqrt(len(filtered_data2)))

DeltaM, unDeltaM = np.mean(filtered_data2), np.std(filtered_data2)/np.sqrt(len(filtered_data2))
print('Number of stars used for magnitude calibration: %s'%len(filtered_data2))

plt.close()
fig = plt.figure(555)
fig.canvas.mpl_connect('button_press_event', onclickclose)
plt.hist(deltaMag, histtype='step', color='r')
plt.hist(filtered_data2, histtype='step', color='b')
plt.legend(['Unfiltered','Filtered'])
plt.xlabel('Delta Mag')
plt.ylabel('Number')
plt.show()

print(band1, band2, band3)
print(band1, K, unK, DeltaM, unDeltaM)

if band1 == 'g': #Make the new file
    f = open(reducedpath+'Corrections.txt', 'w')
    f.write('band K Kun deltaM deltaMun\n')
    f.write('%s %s %s %s %s\n'%(band1, K, unK, DeltaM, unDeltaM))
    f.close()
else:
    f = open(reducedpath+'Corrections.txt', 'a')
    f.write('%s %s %s %s %s\n'%(band1, K, unK, DeltaM, unDeltaM))
    f.close()



