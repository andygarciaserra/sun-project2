##IMPORTANDO PAQUETES Y  VARIABLES GENERALES##
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
from sympy import *
FIGDIR= 'figures/'
RADTODEG = 180/np.pi
C = 4.67e-13
h = 6.63e-34
k = 1.38e-23
c = 3e8
Tm = 6300
#strong_reg = [175, 200, 225, 250]          #Región de Umbra Andy
#strong_reg = [210, 230, 170, 200]          #Región de Umbra Andy NUEVA
#weak_reg = [25, 75, 325, 375]              #Región de sol en calma Andy
weak_reg = [335, 395, 5, 65]               #Región de sol en calma Carlos
strong_reg = [210, 230, 170, 190]          #Región de Umbra Carlos
#strong_reg = [180, 190, 220, 230]          #Región de Umbra Sergio
#weak_reg = [0, 30, 0, 30]                  #Región de sol en calma Sergio



#Creando array dupla de un parámetro de stokes y lamda
def Duple_array(I,lamda,B):
    if(B=='weak'):
        dup_I = np.vstack((lamda,np.mean(I[weak_reg[0]:weak_reg[1],weak_reg[2]:weak_reg[3],:],axis=(0,1)))).T
    elif(B=='strong'):
        dup_I = np.vstack((lamda,np.mean(I[strong_reg[0]:strong_reg[1],strong_reg[2]:strong_reg[3],:],axis=(0,1)))).T
    else:
        print("La región '"+B+"' no es válida, escribe 'weak' para B débil o 'strong' para B fuerte.")
    return dup_I

#Ploteando una imagen fits
def Plot_fits(img,color):
    plt.figure()
    plt.imshow(img,cmap=color)
    plt.colorbar()
    plt.show()

#Plot simple de los espectros
def Simple_plot(I):
    plt.figure()
    plt.plot(I[:,0],I[:,1])
    plt.show()

#Ploteo con índices de vector
def Range_Plot(I):
    plt.figure()
    plt.plot(range(96),I[:,1])
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
    minfit_pos = Find_minmax(min_fit,min_cont)
    minfit_y = Parabola(min_fit,min_cont[minfit_pos])
    return [min_cont[minfit_pos],minfit_y]

#Encontrando máximo dada longitud de onda y rango
def Maximum_Finder(spectra,delta,surr):
    max_pos = np.argmax(spectra[surr[0]:surr[1],1]) + surr[0]
    I_max_surr = spectra[max_pos-delta:max_pos+delta,1]
    x_max_surr = spectra[max_pos-delta:max_pos+delta,0]
    max_fit = np.polyfit(x_max_surr,I_max_surr,2)
    max_cont = np.linspace(spectra[max_pos-delta,0],spectra[max_pos+delta,0],1000)
    maxfit_pos = Find_minmax(max_fit,max_cont)
    maxfit_y = Parabola(max_fit,max_cont[maxfit_pos])
    return [max_cont[maxfit_pos],maxfit_y]

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
        plt.tight_layout()
        plt.savefig(FIGDIR+'Calibration.png', dpi=150)
        plt.show()
    return hinode[:,0]

#Función parábola dado fit y eje x
def Parabola(fit, x):
    return fit[0]*(x**2) + fit[1]*x + fit[2]

#Buscando mínimos dado fit y lamda con la derivada del fit en cada punto
def Find_minmax(fit, xcoord):
    def df(x):
        return 2*fit[0]*x + fit[1]
    dy = np.absolute(df(xcoord))
    minpos = np.argmin(dy)
    return minpos

#Curva (recta) de calibración dados los dos puntos
def Calibration_Curve(p1,p2):
    return np.polyfit(p1,p2,1)

#Calculando y ploteando anchuras Zeeman y Doppler en I+V
def Widenings(I,IplusV,IminV,plot):
    min1 = Minimum_Finder(I,4,[20,30])
    IplusV_cont = np.mean(IplusV[80:87,1])
    min_plus = Minimum_Finder(IplusV,3,[60,70])
    min_min = Minimum_Finder(IminV,3,[70,80])
    hmin = (IplusV_cont + min_plus[1])/2
    pos_izq = np.argmin(IplusV[57:63,1]-hmin)+57
    y_izq_fit = IplusV[pos_izq-2:pos_izq+2,1]
    x_izq_fit = IplusV[pos_izq-2:pos_izq+2,0]
    fit = np.polyfit(y_izq_fit,x_izq_fit,2)
    xizq = Parabola(fit,hmin)
    lam_D = min_plus[0]-xizq
    lam_B = np.absolute(min_min[0]-min_plus[0])/2

    B = lam_B/(C*2.5*(min1[0]**2))

    if(plot=='yes'):
        plt.figure(figsize=(8,6))
        plt.plot(IplusV[:,0],IplusV[:,1])
        plt.plot(IminV[:,0],IminV[:,1])
        plt.scatter(xizq,hmin,color='red',s=30,marker='o',zorder=3)
        plt.scatter(min_plus[0],hmin,color='red',s=30,marker='o',zorder=3)
        plt.hlines(hmin,xmin=xizq,xmax=min_plus[0],color='red', \
                  label=r'$\Delta \lambda_D$ = '+'{0:.2f}'.format(lam_D)+r' $\AA$')
        plt.scatter(min_min[0],min_min[1],color='black',s=30,marker='o',zorder=3)
        plt.scatter(min_plus[0],min_plus[1],color='black',s=30,marker='o',zorder=3)
        plt.hlines(min_plus[1],xmin=min_plus[0],xmax=min_min[0],color='black', \
                  label=r'$\Delta \lambda_B$ = '+'{0:.2f}'.format(lam_B)+r' $\AA$')

        plt.title('B = '+'{0:.2f}'.format(B)+' G')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.savefig(FIGDIR+'Widenings.png', dpi=150)
        plt.show()

    return lam_B, B, min1

#Cálculo de gamma dados Q,U,V
def Gamma(Q,U,V,plot):
    C = 1/(np.sqrt((Q)**2+(U)**2)/V)
    gamma = RADTODEG*np.arccos(np.divide(-1+np.sqrt(4*(C**2)+1),2*C))
    gamma[np.isnan(gamma)] = 0
    hdu = fits.PrimaryHDU(gamma)
    hdu.writeto('gamma.fits')
    gamma = fits.getdata('gamma.fits',ext=0)
    os.remove('gamma.fits')

    if(plot=='yes'):
        Plot_fits(gamma,'turbo')

    return gamma

#Cálculo de phi dados Q,U,V
def Phi(Q,U,V,plot):
    phi = 0.5*RADTODEG*np.arctan2(U,Q)
    phi[np.isnan(phi)] = 0
    phi += 90
    hdu = fits.PrimaryHDU(phi)
    hdu.writeto('phi.fits')
    phi = fits.getdata('phi.fits',ext=0)
    os.remove('phi.fits')
    
    if(plot=='yes'):
        Plot_fits(phi,'turbo')

    return phi

#Cálculo de B_long dados I,V,lam0 y geff
def BLong_Calc(V,I,lam0,g):
    dI = np.divide((I[1:,1]-I[:-1,1]),((I[1:,0]-I[:-1,0])))
    Blong = -(1 / (C * lam0**2 * g)) * \
             ((np.sum(np.multiply(V[1:,1],dI))) / (np.sum(np.multiply(dI,dI))))
    sigma = np.divide(np.std(V[1:,1]),(C * lam0**2 * g * np.sqrt(np.sum(np.multiply(dI,dI)))))
    return Blong, sigma

#Función de Planck con T
def T(I,Icont,lamda,plot):
    T = 1/(1/Tm-(((k*lamda*1e-10)/(h*c))*np.log(I/Icont)))
    
    if(plot=='yes'):
        Plot_fits(T,'turbo')

    return T

#Cálculo de velocidades 2D
def Velocities(I,V,Imin_c,lamda):
    Delta_bin = np.zeros((len(V),len(V[0])))
    v = np.zeros((len(V),len(V[0])))
    Vrange = 300
    for i in range(len(V)):
        for j in range(len(V[0])):
            Ip = np.vstack((lamda,I[i,j,:])).T
            Vp = np.vstack((lamda,V[i,j,:])).T
            max = np.max(Vp[:,1])

            if(max > Vrange):
                delta_lam = Middle_minmax(Vp,Imin_c)
                v[i,j] = (delta_lam*c)/(Imin_c[0])
                Delta_bin[i,j] = 1

            elif(max <= Vrange):
                Imin = Minimum_Finder(Ip,4,[4,50])
                delta_lam = (Imin[0]-Imin_c[0])
                v[i,j] = (delta_lam*c)/(Imin_c[0])
                Delta_bin[i,j] = 0

    return Delta_bin, v

#Cálculo del punto medio entre minimo y máximo
def Middle_minmax(Vp, Imin_c):
    min = Minimum_Finder(Vp,4,[4,50])
    max = Maximum_Finder(Vp,4,[4,50])
    mid = (min[0]+max[0])/2
    delta_lam = mid - Imin_c[0]
    return delta_lam