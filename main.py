##IMPORTING PACKAGES##
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
DATADIR = '/home/andy/Downloads/'
FIGDIR = 'figures/'





#Espectra calibration
    #Importing data from fits
data = fits.getdata(DATADIR+'Stokes_sunspot_HINODE.fits')
    
    #Creating raw intensity image at t=1
if (os.path.exists('new.fits')):
    os.remove('new.fits')
hdu = fits.PrimaryHDU(data[0,:,:,50])
hdu.writeto('new.fits')
I_img = fits.getdata('new.fits',ext=0)
plt.figure()
plt.imshow(I_img, cmap='viridis')
plt.colorbar()
plt.tight_layout()
plt.savefig(FIGDIR+'Full_Image.png', dpi=150)
plt.show()
    
    #Plotting spectra for calm region
I = data[0,:,:,:]        #Calm region at (350,50)
I_calm = I[25:75,325:375,:]
I_dark = I[175:200,225:250,:]
mean_calm = np.mean(I_calm, axis=(0,1))
mean_dark = np.mean(I_dark, axis=(0,1))
plt.figure()
plt.plot(range(96),mean_calm)
plt.plot(range(96),mean_dark)
plt.tight_layout()
plt.savefig(FIGDIR+'Calm_vs_Dark_t50.png', dpi=150)
plt.show()

    #Atlas image calibrating
atlas = fits.getdata('data/atlas_6301_6302.fits')
# x_min1 = np.min(atlas[740:780,0])
y_min1 = np.min(atlas[740:780,1])
# x_min2 = np.min(atlas[1225:1275,0])
y_min2 = np.min(atlas[1225:1275,1])
# print(x_min1,y_min1)
plt.figure()
plt.hlines(y_min1,xmin = np.min(atlas[:,0]),xmax = np.max(atlas[:,0]),color='red',ls='--')
plt.hlines(y_min2,xmin = np.min(atlas[:,0]),xmax = np.max(atlas[:,0]),color='orange',ls='--')
plt.plot(atlas[:,0],atlas[:,1])
plt.show()

