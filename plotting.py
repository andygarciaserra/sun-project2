##IMPORTANDO PAQUETES Y VARIABLES GENERALES##
from func import *
import warnings
warnings.filterwarnings("ignore")
FIGDIR = 'figures/'
DATADIR = '/home/andy/Downloads/prac2solar/data/'

data = Import_data(DATADIR+'Stokes_sunspot_HINODE.fits')
lamda = Img_calibration(DATADIR, 'no')
Iw = Duple_array(data[0,:,:,:],lamda,'weak')
Is = Duple_array(data[0,:,:,:],lamda,'strong')
factor = Iw[0,1]/Is[0,1]
print(factor)
plt.figure()
plt.plot(Iw[:,0],Iw[:,1], label = 'Sol en calma')
plt.plot(Is[:,0],factor*Is[:,1], label = 'Umbra')
plt.xlabel('longitud de onda ['+r'$\AA$'+']')
plt.ylabel('I [cuentas]')
plt.legend()
plt.tight_layout()
plt.savefig(FIGDIR+'Iw_vs_Is.png', dpi=150)
plt.show()
