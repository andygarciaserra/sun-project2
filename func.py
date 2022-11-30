##IMPORTANDO PAQUETES##
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from sympy import *





#Plot simple de los espectros
def Simple_plot(I):
    plt.figure()
    plt.plot(I[:,0],I[:,1])
    plt.show()

#Importando data de un fits
def Import_data(dir):
    data = fits.getdata(dir)
    return data

#Visualizando mancha en t=50
def print_img(figdir):
    data = Import_data(figdir)
    hdu = fits.PrimaryHDU(data[0,:,:,50])
    hdu.writeto('new.fits')
    I_img = fits.getdata('new.fits',ext=0)
    os.remove('new.fits')
    plt.figure()
    plt.imshow(I_img, cmap='viridis')
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(figdir+'Full_Image.png', dpi=150)
    plt.show()

#Calculando espectros medios de las intensidades
def Mean_Spectra_regions(FIGDIR,I_calm,I_dark,plot):
    mean_calm = np.mean(I_calm, axis=(0,1))
    mean_dark = np.mean(I_dark, axis=(0,1))
    if (plot=='yes'):
        plt.figure()
        plt.plot(range(96),mean_calm)
        plt.plot(range(96),mean_dark)
        plt.tight_layout()
        plt.savefig(FIGDIR+'Calm_vs_Dark_t50.png', dpi=150)
        plt.show()
    return(mean_calm,mean_dark)

#Encontrando dos mínimos dada longitud de onda y rango
def Double_Minimum_Finder(spectra,delta,surr_min1,surr_min2):
    xmin1,ymin1 = Minimum_Finder(spectra,delta,surr_min1)
    xmin2,ymin2 = Minimum_Finder(spectra,delta,surr_min2)
    return [xmin1,xmin2,ymin1,ymin2]

#Encontrando mínimo dada longitud de onda y rango
def Minimum_Finder(spectra,delta,surr):
    min_pos = np.argmin(spectra[surr[0]:surr[1],1]) + surr[0]
    I_min_surr = spectra[min_pos-delta:min_pos+delta,1]
    x_min_surr = spectra[min_pos-delta:min_pos+delta,0]
    min_fit = np.polyfit(x_min_surr,I_min_surr,2)
    min_cont = np.linspace(spectra[min_pos-delta,0],spectra[min_pos+delta,0],1000)
    minfit_pos = Find_minimum(min_fit,min_cont)
    minfit_y = min_fit[0]*min_cont[minfit_pos]**2 + min_fit[1]*min_cont[minfit_pos]+min_fit[2]
    return [min_cont[minfit_pos],minfit_y]


# Calibracion de los datos HINODE usando el espectro de referncia de atlas
def Img_calibration(DATADIR, plot):
    
    #Importando datos
    atlas = Import_data(DATADIR+'atlas_6301_6302.fits')
    data = Import_data(DATADIR+'Stokes_sunspot_HINODE.fits')
    hinode = np.vstack((range(96),np.mean(data[0,25:75,325:375,:],axis=(0,1)))).T

    #Creando vectores de alrededores de atlas para calcular mínimos:       
    xcoords_atlas = Double_Minimum_Finder(atlas,15,[740,780],[1225,1275])
    xcoords_hinode = Double_Minimum_Finder(hinode,8,[15,35],[60,80])
    min1_factor = hinode[0,1]/atlas[0,1]
    min2_factor = hinode[0,1]/atlas[0,1]
    vfactor = np.mean([min1_factor,min2_factor])

    #Calibración final
    curve = Calibration_Curve([xcoords_hinode[0], xcoords_hinode[1]],
                              [xcoords_atlas[0], xcoords_atlas[1]])
    atlas[:,1] = vfactor*atlas[:,1]
    hinode[:,0] = hinode[:,0]*curve[0]+curve[1]

    #Plots
    if (plot=='yes'):
        plt.figure()
        plt.plot(atlas[:,0],atlas[:,1],label='Atlas data')
        plt.plot(hinode[:,0],hinode[:,1],label='HINODE calibration')
        plt.legend()
        plt.show()
    return hinode[:,0]

    #Función parábola dado fit y eje x
def Parabola(fit, x):
    return fit[0]*(x**2) + fit[1]*x + fit[2]

    #Buscando mínimos dado fit y lamda con la derivada del fit en cada punto
def Find_minimum(fit, xcoord):
    def df(x):
        return 2*fit[0]*x + fit[1]
    dy = np.absolute(df(xcoord))
    minpos = np.argmin(dy)
    return minpos

    #Curva (recta) de calibración dados los dos puntos
def Calibration_Curve(p1,p2):
    return np.polyfit(p1,p2,1)
