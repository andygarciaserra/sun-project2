##IMPORTANDO PAQUETES Y VARIABLES GENERALES##
from func import *
import warnings
warnings.filterwarnings("ignore")
FIGDIR = 'figures/'
DATADIR = '/home/andy/Downloads/prac2solar/data/'


    #2# CALIBRACION:
#print_img(DATADIR+'Stokes_sunspot_HINODE.fits')
lamda = Img_calibration(DATADIR, 'no')


    #3# MEDIDA DEL CAMPO MAGNÉTICO:
        #3.1# Campo fuerte:
        #ploteando I para ver desdoblamiento
data = Import_data(DATADIR+'Stokes_sunspot_HINODE.fits')

        #calculando anchuras Zeeman, Doppler e intensidad de B en I+V
I = Duple_array(data[0,:,:,:],lamda,'weak')
IplusV = Duple_array(np.add(data[0,:,:,:],data[3,:,:,:]),lamda,'strong')
IminV = Duple_array(data[0,:,:,:]-data[3,:,:,:],lamda,'strong')
lam_B, B, min1 = Widenings(I,IplusV,IminV,'yes')

        #importando 2D de Q,U,V
Q = data[1,:,:,67]
U = data[2,:,:,67]
V = data[3,:,:,67]

        #calculando phi y ploteando
#phi = Phi(Q,U,V,'yes')

        #calculando gamma y ploteando
#gamma = Gamma(Q,U,V,'yes')

        #3.2#Campo débil:
        #calculando B_long
I = Duple_array(data[0,:,:,:],lamda,'weak')
V = Duple_array(data[3,:,:,:],lamda,'weak')
lam0, geff = min1[0], 1.667
B_long, sigma = BLong_Calc(V,I,lam0,geff)
print('B= '+'{0:.2f}'.format(B_long)+' +/- '+'{0:.2f}'.format(sigma)+' Gauss')

    #4# MAPA DE TEMPERATURA Y CAMPO DE VELOCIDADES
        #calculo de Temperatura
Icont = np.mean(I[40:60,1])
I = data[0,:,:,67]
T = T(I,Icont,lamda[67],'yes')

        #velocidades
I = Duple_array(np.add(data[0,:,:,:],data[3,:,:,:]),lamda,'weak')
Imin_c = Minimum_Finder(I,2,[0,96])
I = data[0,:,:,:]
V = data[3,:,:,:]
Delta_lam, v = Velocities(I,V,Imin_c,lamda)

hdu = fits.PrimaryHDU(Delta_lam)
hdu.writeto('delta.fits')
delta = fits.getdata('delta.fits',ext=0)
os.remove('delta.fits')
Plot_fits(delta,'turbo')

v = v*1e-3
plt.figure()
plt.imshow(v,cmap='bwr',vmin=-2,vmax=1.5)
plt.xlabel('coordenada X [pix]')
plt.ylabel('coordenada Y [pix]')
plt.colorbar(label='v [km/s]')
plt.tight_layout()
plt.savefig(FIGDIR+'2D_v".png', dpi=150)
plt.show()