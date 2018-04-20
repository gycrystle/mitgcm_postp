# python

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


dyU = 0.0*np.ones((250,250));
dxV = 0.0*np.ones((250,250));

mld  = 0.0*np.ones((250,250));
mld_mid = 0.0*np.ones((250,250));
mld_below = 0.0*np.ones((250,250));
mld_ind  = 0.0*np.ones((250,250));
w_mld  = 0.0*np.ones((250,250));
w_mld_mid = 0.0*np.ones((250,250));
w_mld_below = 0.0*np.ones((250,250));
         
zeta = 0.0*np.ones((250,250));
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

itrs = [720, 1440, 2160, 2880, 3600] #, 2160, 2880, 3600, 4320, 5040, 5760, 6480, 7200]
for it in itrs:
    U = mit.rdmds('U',it)
    V = mit.rdmds('V',it)
    W = mit.rdmds('W',it)
    Kzu = mit.rdmds('GGL90viscArU', it)
    Kzv = mit.rdmds('GGL90viscArV', it)

    taux = mit.rdmds('diagTAUX',it)
    tauy = mit.rdmds('diagTAUY',it)
    
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
            dyU[j,i] = (U[1,j,i]-U[1,j-1,i])/dYC[j,i];
            dxV[j,i] = (V[1,j,i]-V[1,j,i-1])/dXC[j,i];
            zeta[j,i] = dxV[j,i]-dyU[j,i];
#
    for i in range(1,249):
        for j in range(1,249):
            zeta_intx[j,i] =(zeta[j,i]+zeta[j+1,i])/2
            zeta_inty[j,i] =(zeta[j,i]+zeta[j,i+1])/2
#
            termx[j,i] =taux[j,i]/(1e-4+zeta_intx[j,i])
            termy[j,i] =tauy[j,i]/(1e-4+zeta_inty[j,i])
#
    for i in range(1,248):
        for j in range(1,248):
            dytermx[(j+1),i] = (termx[j+1,i]-termx[j,i])/dYC[j+1,i];
            dxtermy[j,(i+1)] = (termy[j,i+1]-termy[j,i])/dXC[j,i+1];
            W_tot[(j+1),(i+1)] = (dxtermy[j+1,i+1]-dytermx[j+1,i+1])/rho0;


    #plt.figure()
    fig = plt.figure(figsize=(10,16))
    
    
    ax1 = fig.add_subplot(4, 2, 1)
    ax1.set_aspect(1)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf(XC*1e-3,YC*1e-3,mld_mid, mld_range,cmap=cm.Blues_r)
    plt.colorbar(label="depth (m)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("Mid Mixed Layer")
    #plt.savefig("mld"+ str(i) + ".png")

    ax7 = fig.add_subplot(4, 2, 3)
    ax7.set_aspect(1)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf(XC*1e-3,YC*1e-3,mld, mld_range,cmap=cm.Blues_r)
    plt.colorbar(label="depth (m)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("Mixed Layer Depth")

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

    
    plt.tight_layout(pad=1)
    plt.savefig("../figures_relwind/w_mld2_ekmanpump_"+ str(it) + ".png")


#plt.contour(xi, yi, zi, v, linewidths=0.5, colors='k')
#plt.contourf(xi, yi, zi, v, cmap=plt.cm.jet)
#x = plt.colorbar(ticks=v)


