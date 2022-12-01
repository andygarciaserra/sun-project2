from func import *
FIGDIR = 'figures/'
DATADIR = '/home/andy/Downloads/prac2solar/data/'
RADTODEG = 180/np.pi


# Las regiones en uso:
#       - Sol en calma: 25:75,325:375
#       - Umbra: 175:200,225:250



    #2# CALIBRACION:
#print_img(DATADIR+'Stokes_sunspot_HINODE.fits')
lamda = Img_calibration(DATADIR, plot='no')



    #3# MEDIDA DEL CAMPO MAGNÉTICO:
        #3.1# Campo fuerte:
        #ploteando I para ver desdoblamiento
data = Import_data(DATADIR+'Stokes_sunspot_HINODE.fits')
#I = np.vstack((lamda,np.mean(data[0,175:200,225:250,:],axis=(0,1)))).T
#Simple_plot(I)

        #calculando anchuras Zeeman y Doppler en I-V
I_V = np.vstack((lamda,np.add(np.mean(data[0,175:200,225:250,:],axis=(0,1)),
                 np.mean(data[3,175:200,225:250,:],axis=(0,1))).T)).T

I_V_cont = np.mean(I_V[33:47,1])
min1 = Minimum_Finder(I_V,2,[10,30])
hmin = (I_V_cont - min1[1])/2 + min1[1]
min = np.argmin(I_V[:,1])
izq = np.argmin(np.absolute(I_V[0:min,1]-hmin))
fit_izq = np.polyfit(I_V[izq-2:izq+1,0],I_V[izq-2:izq+1,1],2)
x_izq = np.linspace(I_V[izq-2,0],I_V[izq+1,0],10000)
y_izq = Parabola(fit_izq,x_izq)
pos_izq = np.argmin(np.absolute(y_izq-hmin))
min2 = Minimum_Finder(I_V,2,[62,68])
min3 = Minimum_Finder(I_V,2,[69,72])
lam_D = min1[0]-x_izq[pos_izq]
lam_B = min3[0]-min2[0]

        #intensidad de B
B = lam_B/((4.67e-13)*2.5*(min1[0]**2))

        #plotting
plt.figure(figsize=(8,6))
plt.plot(I_V[:,0],I_V[:,1])
plt.plot(x_izq,y_izq,zorder=3)
plt.scatter(x_izq[pos_izq],y_izq[pos_izq],color='red',s=30,marker='o',zorder=3)
plt.scatter(I_V[izq-2,0],I_V[izq-2,1],color='orange',s=30,marker='+',zorder=2)
plt.scatter(I_V[izq+1,0],I_V[izq+1,1],color='orange',s=30,marker='+',zorder=2)
plt.scatter(min3[0],min2[1],color='black',s=15,marker='o',zorder=3)
plt.scatter(min1[0],hmin,color='orange',s=30,marker='+',
            label='{0:.2f}'.format(hmin),zorder=2)
plt.scatter(min1[0],min1[1],color='blue',s=20,
            label='{0:.2f}'.format(min1[0])+', '+'{0:.2f}'.format(min1[1]),
            zorder=3)
plt.scatter(min2[0],min2[1],color='orange',s=20,
            label='{0:.2f}'.format(min2[0])+', '+'{0:.2f}'.format(min2[1]),
            zorder=3)
plt.scatter(min3[0],min3[1],color='purple',s=20,
            label='{0:.2f}'.format(min3[0])+', '+'{0:.2f}'.format(min3[1]),
            zorder=3)
plt.vlines(min3[0], ymin=min2[1], ymax=min3[1],ls='--',color='k', alpha=0.4)
plt.hlines(min2[1], xmin=min2[0], xmax=min3[0],ls='-',color='k',
           label=r'$\Delta \lambda_B$ = '+'{0:.2f}'.format(lam_B)+r' $\AA$')
plt.hlines(hmin, xmin=x_izq[pos_izq], xmax=min1[0],ls='-',color='k',
           label=r'$\Delta \lambda_D$ = '+'{0:.2f}'.format(lam_D)+r' $\AA$')
plt.hlines(I_V_cont,xmin=np.min(I_V[:,0]),xmax=np.max(I_V[:,0]),color='r',
           ls='--',label='{0:.2f}'.format(I_V_cont))
plt.title('B = '+'{0:.2f}'.format(B)+' G')
plt.legend(loc='best')
plt.show()

        #importando 2D de Q,U,V
Q = data[1,:,:,67]
U = data[2,:,:,67]
V = data[3,:,:,67]

        #calculando phi y ploteando
phi = 0.5*RADTODEG*np.arctan2(U,Q)
phi[np.isnan(phi)] = 0
hdu = fits.PrimaryHDU(phi)
hdu.writeto('phi.fits')
phi = fits.getdata('phi.fits',ext=0)
os.remove('phi.fits')
plt.figure()
plt.imshow(phi,cmap='bwr')
plt.colorbar()
plt.show()

        #calculando gamma y ploteando
C = 1/(np.sqrt((Q)**2+(U)**2)/V)
gamma = RADTODEG*np.arccos(np.divide(-1+np.sqrt(4*(C**2)+1),2*C))
gamma[np.isnan(gamma)] = 0
hdu = fits.PrimaryHDU(gamma)
hdu.writeto('gamma.fits')
gamma = fits.getdata('gamma.fits',ext=0)
os.remove('gamma.fits')
plt.figure()
plt.imshow(gamma,cmap='bwr')
plt.colorbar()
plt.show()





