from func import *
FIGDIR = 'figures/'
DATADIR = '/home/andy/Downloads/prac2solar/data/'

# We are using the regions:
#       - Calm Sun: 25:75,325:375
#       - Dark spot: 175:200,225:250

    #2# Calibration:
#print_img(DATADIR+'Stokes_sunspot_HINODE.fits')
lamda = Img_calibration(DATADIR, plot='no')

    #3.1.1# Magnetic field:
data = Import_data(DATADIR+'Stokes_sunspot_HINODE.fits')
#I_V = np.vstack((lamda,np.add(np.mean(data[0,175:200,225:250,:],axis=(0,1)),
#                 np.mean(data[3,175:200,225:250,:],axis=(0,1))).T)).T
#print(I)
#I_V_cont = np.mean(I_V[33:47,1])
#min1 = Minimum_Finder(I_V,2,[10,30])
#hmin = (I_V_cont - min1[1])/2 
#
#    #Plotting
#plt.figure()
#plt.plot(I_V[:,0],I_V[:,1])
#plt.scatter(min1[0],hmin+min1[1],color='orange',s=30,marker='+',
#            label='{0:.2f}'.format(hmin),zorder=2)
#plt.scatter(min1[0],min1[1],color='orange',s=15,
#            label='{0:.2f}'.format(min1[0])+', '+'{0:.2f}'.format(min1[1]),
#            zorder=1)
#plt.hlines(I_V_cont,xmin=np.min(I_V[:,0]),xmax=np.max(I_V[:,0]),color='r',
#           ls='--',label='{0:.2f}'.format(I_V_cont))
#plt.legend()
#plt.show()

    #3.1.3# Acimuth
Q = data[1,:,:,:]
U = data[2,:,:,:]
phi = 0.5*np.arctan2(U,Q,)