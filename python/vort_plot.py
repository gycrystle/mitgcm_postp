"""
Python script to read and plot vertical velocity (W) output from MITgcm
- Define Mixed layer depth with deltaT>0.5°C from the surface
- Plot W at Mixed Layer depth
- Plot Mid MLD and it's respective W
- Plot W below mixed layer depth
- Plot the wind stress
- Compute Theoretical Ekman pumping
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm

plt.ion()

rho0 = 1028.0
XC = mit.rdmds('XC')
dXC = mit.rdmds('DXC')
YC = mit.rdmds('YC')
dYC = mit.rdmds('DYC')
RC = mit.rdmds('RC')

dyU = 0.0*np.ones((40,250,250));
dxV = 0.0*np.ones((40,250,250));

mld  = 0.0*np.ones((250,250));
mld_mid = 0.0*np.ones((250,250));
mld_below = 0.0*np.ones((250,250));
mld_ind  = 0.0*np.ones((250,250)); #Mixed layer depth indices
w_mld  = 0.0*np.ones((250,250));
w_mld_mid = 0.0*np.ones((250,250));
w_mld_below = 0.0*np.ones((250,250));
         
zeta = 0.0*np.ones((40,250,250));
zeta_intx = 0.0*np.ones((250,250));
zeta_inty = 0.0*np.ones((250,250));
#
termx =0.0*np.ones((250,250));
termy =0.0*np.ones((250,250));
#
dytermx = 0.0*np.ones((250,250));
dxtermy = 0.0*np.ones((250,250));
W_tot = 0.0*np.ones((250,250));

mld_range = np.linspace(-275.0, -25, 26, endpoint=True)
w_range = np.linspace(-0.5, 0.5, 101, endpoint=True)

itrs = [2160, 2880, 3600] #, 2160, 2880, 3600, 4320, 5040, 5760, 6480, 7200]
for it in itrs:
    U = mit.rdmds('U',it)
    V = mit.rdmds('V',it)
    W = mit.rdmds('W',it)
#    taux = mit.rdmds('diagTAUX',it)
#    tauy = mit.rdmds('diagTAUY',it)
    
    Temp = mit.rdmds('T',it)

    temp_surf= Temp[1,:,:]
    
    Tdiff=abs(Temp-temp_surf)

    for x in range(0, 250):
        for y in range(0, 250):
            mld_ind[y,x]=np.where(Tdiff[:,y,x] > 0.5)[0][0]
            mld[y,x]=RC[int(mld_ind[y,x])]
            mld_mid[y,x]=RC[int(round(mld_ind[y,x]/2))]
            mld_below[y,x]=RC[int(mld_ind[y,x]+9)]
            w_mld[y,x]=W[int(mld_ind[y,x]),y,x]
            w_mld_mid[y,x]=W[int(round(mld_ind[y,x]/2)),y,x]
            w_mld_below[y,x]=W[int(mld_ind[y,x])+9,y,x]

    for i in range(1,250):
        for j in range(1,250):
            #for k in range(0,40):
            dyU[:,j,i] = (U[:,j,i]-U[:,j-1,i])/dYC[j,i];
            dxV[:,j,i] = (V[:,j,i]-V[:,j,i-1])/dXC[j,i];
            zeta[:,j,i] = dxV[:,j,i]-dyU[:,j,i];
#
#Plot all figures
    #plt.figure()
    fig = plt.figure(figsize=(12,4))
    
    
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf(XC*1e-3,YC*1e-3,zeta[1,:,:]*1e5, 100,cmap=cm.inferno)
    plt.colorbar(label=r'$\zeta$')
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("Vorticity")
    #plt.savefig("mld"+ str(i) + ".png")

    ax7 = fig.add_subplot(1, 2, 2)
    ax7.set_aspect(0.1)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf(XC[125,:]*1e-3,RC.squeeze(),zeta[:,125,:]*1e5, 100,cmap=cm.hsv)
    plt.colorbar(label="depth (m)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("Mixed Layer Depth")
    """
    ax8 = fig.add_subplot(4, 2, 5)
    ax8.set_aspect(1)
    #plt.subplot(121, aspect='equal', adjustable='box-forced)
    plt.contourf(XC*1e-3,YC*1e-3,mld_below, mld_range,cmap=cm.Blues_r)
    plt.colorbar(label="depth (m)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("Below Mixed Layer")


    #plt.figure()
    ax4 = fig.add_subplot(4, 2, 2)
    ax4.set_aspect(1)
    #plt.subplot(122, aspect='equal', adjustable='box-forced')
    plt.contourf(XC*1e-3,YC*1e-3,w_mld_mid*1e3, w_range,cmap=cm.seismic)
    plt.colorbar(label="W (mm/s)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("W in the middle of mixed layer")
    
    #plt.figure()
    ax5 = fig.add_subplot(4, 2, 4)
    ax5.set_aspect(1)
    #plt.subplot(122, aspect='equal', adjustable='box-forced')
    plt.contourf(XC*1e-3,YC*1e-3,w_mld*1e3, w_range,cmap=cm.seismic)
    plt.colorbar(label="W (mm/s)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("W at the base of mixed layer")

    #plt.figure()
    ax6 = fig.add_subplot(4, 2, 6)
    ax6.set_aspect(1)
    #plt.subplot(122, aspect='equal', adjustable='box-forced')
    plt.contourf(XC*1e-3,YC*1e-3,w_mld_below*1e3, w_range,cmap=cm.seismic)
    plt.colorbar(label="W (mm/s)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("W below mixed layer")

    
    ax3 = fig.add_subplot(4, 2, 7)
    ax3.set_aspect(1)
    #plt.subplot(122, aspect='equal', adjustable='box-forced')
    plt.contourf(XC*1e-3,YC*1e-3,taux)
    plt.colorbar(label="TAUx (N/m^2)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("Wind stress")

    ax2 = fig.add_subplot(4, 2, 8)
    ax2.set_aspect(1)
    #plt.subplot(122, aspect='equal', adjustable='box-forced')
    plt.contourf(XC*1e-3,YC*1e-3,W_tot*1e3, w_range,cmap=cm.seismic)
    plt.colorbar(label="W_tot (mm/s)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("Ekman pumping")
    """
    
    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/vort_000"+ str(it) + ".png")
    elif it < 1000:
        plt.savefig("./figures/vort_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/vort_"+ str(it) + ".png")
#    plt.close()


#plt.contour(xi, yi, zi, v, linewidths=0.5, colors='k')
#plt.contourf(xi, yi, zi, v, cmap=plt.cm.jet)
#x = plt.colorbar(ticks=v)


