##IMPORTING PACKAGES##
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
DATADIR = '/home/andy/Downloads/'



#Importing data from fits
data = fits.getdata(DATADIR+'Stokes_sunspot_HINODE.fits')
I = data[0,:,:,:]
Q = data[1,:,:,:]
U = data[2,:,:,:]
V = data[3,:,:,:]


#Espectra calibration
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
plt.show()
    #Plotting spectra for calm region
I = data[0,350,50,:]        #Calm region at (350,50)
x_axis = range(len(I))
mean = np.mean(I)
plt.figure()
plt.plot(x_axis,I)
plt.hlines(mean, x_axis[0], x_axis[len(x_axis)-1],'r',ls='--')
plt.tight_layout()
plt.show()

# EN VEZ DE TOMAR UN PUNTO PARA CALCULAR EL VALOR MEDIO SELECCIONAR UNA REGIÓN DE 10-15 PÍXELES DE ANCHO
# Y HACER EL VALOR MEDIO DE LA INTENSIDAD EN ESA REGIÓN PARA CADA LONGITUD DE ONDA, ASÍ NOS QUEDA UN ESPECTRO
# MEDIADO EN INTENSTIDAD POR UNIDAD DE LONGITUD DE ONDA.

