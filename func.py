##IMPORTING PACKAGES##
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from sympy import *





#Espectra calibration

    #Importing data from fits
def Import_data(DATADIR):
    data = fits.getdata(DATADIR+'Stokes_sunspot_HINODE.fits')
    return data

    #Creating raw intensity image at t=50
def print_img(figdir,data):
    if (os.path.exists('new.fits')):
        os.remove('new.fits')
    hdu = fits.PrimaryHDU(data[0,:,:,50])
    hdu.writeto('new.fits')
    I_img = fits.getdata('new.fits',ext=0)
    plt.figure()
    plt.imshow(I_img, cmap='viridis')
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(figdir+'Full_Image.png', dpi=150)
    plt.show()
        
        #Plotting spectra for calm and dark regions
def Import_regions(data):
    I = data[0,:,:,:]           #Calm region at (350,50)
    I_calm = I[25:75,325:375,:]
    I_dark = I[175:200,225:250,:]
    return(I_calm,I_dark)

def Mean_Spectra_regions(FIGDIR,I_calm,I_dark):
    mean_calm = np.mean(I_calm, axis=(0,1))
    mean_dark = np.mean(I_dark, axis=(0,1))
    plt.figure()
    plt.plot(range(96),mean_calm)
    plt.plot(range(96),mean_dark)
    plt.tight_layout()
    plt.savefig(FIGDIR+'Calm_vs_Dark_t50.png', dpi=150)
    plt.show()
    
        #Atlas image calibrating
def Img_calibration(DATADIR):
    dlambda = 30
    atlas = fits.getdata(DATADIR+'atlas_6301_6302.fits')
    min1_pos = np.argmin(atlas[740:780,1]) + 740
    min2_pos = np.argmin(atlas[1225:1275,1]) + 1225
    I_min1_surr = atlas[min1_pos-dlambda:min1_pos+dlambda,1]
    x_min1_surr = atlas[min1_pos-dlambda:min1_pos+dlambda,0]
    I_min2_surr = atlas[min2_pos-dlambda:min2_pos+dlambda,1]
    x_min2_surr = atlas[min2_pos-dlambda:min2_pos+dlambda,0]

    #Fitting to find minimums
    min1_fit = np.polyfit(x_min1_surr,I_min1_surr,2)
    min2_fit = np.polyfit(x_min2_surr,I_min2_surr,2)
    min1fit_pos = Find_minimum(min1_fit,x_min1_surr)
    print(x_min1_surr[min1fit_pos])
    min2fit_pos = Find_minimum(min2_fit,x_min2_surr)
    print(x_min2_surr[min2fit_pos])

    # Plotting
    xmin1_fit = np.linspace(atlas[min1_pos-dlambda,0],atlas[min1_pos+dlambda,0],60)
    ymin1_fit = min1_fit[0]*(xmin1_fit)**2 + min1_fit[1]*(xmin1_fit) + min1_fit[2]
    xmin2_fit = np.linspace(atlas[min2_pos-dlambda,0],atlas[min2_pos+dlambda,0],60)
    ymin2_fit = min2_fit[0]*(xmin2_fit)**2 + min2_fit[1]*(xmin2_fit) + min2_fit[2]
    plt.figure(x_min1_surr[min1fit_pos], ymin=0, ymax = 10000)
    plt.vlines(x_min2_surr[min2fit_pos], ymin=0, ymax = 10000)
    plt.scatter(atlas[min1_pos,0],atlas[min1_pos,1],c='r',s=10)
    plt.scatter(atlas[min2_pos,0],atlas[min2_pos,1],c='r',s=10)
    plt.plot(atlas[:,0],atlas[:,1])
    plt.plot(xmin1_fit, ymin1_fit)
    plt.plot(xmin2_fit, ymin2_fit)
    plt.show()

def Find_minimum(fit,xcoord):
    def df(x):
        return 2*fit[0]*x + fit[1]
    dy = np.absolute(df(xcoord))
    minpos = np.argmin(dy)
    return minpos