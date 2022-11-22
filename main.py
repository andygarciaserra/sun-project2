from func import *
FIGDIR = 'figures/'
DATADIR = '/home/andy/Downloads/prac2solar/data/'


# data = Import_data(DATADIR+'Stokes_sunspot_HINODE.fits')
# print_img(FIGDIR,data)
# I_calm, I_dark = Import_regions(data)
# Mean_Spectra_regions(DATADIR,I_calm, I_dark)
Img_calibration(DATADIR, plot='yes')