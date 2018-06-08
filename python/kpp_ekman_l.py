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
day_s = 2
day_e = 3

dumpfreq = 86400
timestep = 120
nstart = int(day_s*86400/dumpfreq) # integer of dumpfrec, either start from zero or 
nend = int(day_e*86400/dumpfreq) # no of time step 
itrs = dumpfreq/timestep*np.arange(nstart,nend)

rho0 = 999.8
f0 = 1e-4

XC = mit.rdmds('XC')
dXC = mit.rdmds('DXC')
YC = mit.rdmds('YC')
dYC = mit.rdmds('DYC')
RC = mit.rdmds('RC')

nr = RC.size
nx = XC[0,:].size
ny = YC[:,0].size


dyU = np.zeros((ny,nx));
dxV = np.zeros((ny,nx));
mld_ind = np.zeros((ny,nx));
mld = np.zeros((ny,nx));

"""
mld  = 0.0*np.ones((250,250));
mld_mid = 0.0*np.ones((250,250));
mld_below = 0.0*np.ones((250,250));
mld_ind  = 0.0*np.ones((250,250)); #Mixed layer depth indices
w_mld  = 0.0*np.ones((250,250));
w_mld_mid = 0.0*np.ones((250,250));
w_mld_below = 0.0*np.ones((250,250));
"""
         
zeta = np.zeros((ny,nx));
zeta_intx = np.zeros((ny,nx));
zeta_inty = np.zeros((ny,nx));
#
termx = np.zeros((ny,nx));
termy = np.zeros((ny,nx));
#
dytermx = np.zeros((ny,nx));
dxtermy = np.zeros((ny,nx));
W_tot = np.zeros((ny,nx));

zcont = np.linspace(-340,-20,17) # for contour
mld_range = np.linspace(-275.0, -25, 26, endpoint=True)
w_range = np.linspace(-0.15, 0.15, 101, endpoint=True)

#itrs = [720, 1440, 2160, 2880, 3600, 4320, 5040, 5760, 6480, 7200]
for it in itrs:
    U = mit.rdmds('U',it)
    V = mit.rdmds('V',it)
    W = mit.rdmds('W',it)
    taux = mit.rdmds('diagTAUX',it)
    tauy = mit.rdmds('diagTAUY',it)
    nuz = mit.rdmds('KPPviscAz',it)
#    nuz_atv = mit.rdmds('GGL90viscArV',it)
    ek_depth = np.sqrt(2.0*nuz[1,:,:]/f0)
#    ekman_layer_v = np.sqrt(2.0*nuz_atv/f0)
#    
    Temp = mit.rdmds('T',it)
#
    temp_surf= Temp[0,:,:]
#    
    Tdiff=abs(Temp-temp_surf)
#
#
    for x in range(0, nx):
        for y in range(0, ny):
            mld_ind[y,x]=np.where(Tdiff[:,y,x] > 0.5)[0][0]
            mld[y,x]=RC[int(mld_ind[y,x])]
#            mld_mid[y,x]=RC[int(round(mld_ind[y,x]/2))]
#            mld_below[y,x]=RC[int(mld_ind[y,x]+9)]
#            w_mld[y,x]=W[int(mld_ind[y,x]),y,x]
#            w_mld_mid[y,x]=W[int(round(mld_ind[y,x]/2)),y,x]
#            w_mld_below[y,x]=W[int(mld_ind[y,x])+9,y,x]
#
    for i in range(1,nx):
        for j in range(1,ny):
            dyU[j,i] = (U[1,j,i]-U[1,j-1,i])/dYC[j,i];
            dxV[j,i] = (V[1,j,i]-V[1,j,i-1])/dXC[j,i];
            zeta[j,i] = dxV[j,i]-dyU[j,i];
#
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            zeta_intx[j,i] =(zeta[j,i]+zeta[j+1,i])/2
            zeta_inty[j,i] =(zeta[j,i]+zeta[j,i+1])/2
#
            termx[j,i] =taux[j,i]/(1e-4+zeta_intx[j,i])
            termy[j,i] =tauy[j,i]/(1e-4+zeta_inty[j,i])
#
    for i in range(1,nx-2):
        for j in range(1,ny-2):
            dytermx[(j+1),i] = (termx[j+1,i]-termx[j,i])/dYC[j+1,i];
            dxtermy[j,(i+1)] = (termy[j,i+1]-termy[j,i])/dXC[j,i+1];
            W_tot[(j+1),(i+1)] = (dxtermy[j+1,i+1]-dytermx[j+1,i+1])/rho0;
#
#Plot all figures
    #plt.figure()
    fig = plt.figure(figsize=(5,8))
    #plot wind stress
    ax1 = plt.subplot(211)
    ax1.set_aspect(1)
    #
    plt.contourf(XC*1e-3,YC*1e-3,taux)
    plt.colorbar(label=r'$\tau_x \ [N/m^2]$')
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(r'Wind stress $\tau_x$')

    
    #eddy visc
    ax2 = plt.subplot(212)
    ax2.set_aspect(1)
    #
    plt.contourf(XC*1e-3,YC*1e-3,nuz[1,:,:], 100)
    plt.colorbar(label=r'$\nu_z \ [m^2/s]$')
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("Eddy viscosity @ surface")
    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/windstress_visc_0000"+ str(it) + ".png")
    elif it < 100:
        plt.savefig("./figures/windstress_visc_000"+ str(it) + ".png")
    elif it < 1000:
        plt.savefig("./figures/windstress_visc_00"+ str(it) + ".png")
    elif it < 10000:
        plt.savefig("./figures/windstress_visc_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/windstress_visc_"+ str(it) + ".png")
#    plt.close()

    #plot ekman layer depth vs mld
    fig = plt.figure(figsize=(5,8))
    ax1 = plt.subplot(211)
    ax1.set_aspect(1)
#
    plt.contourf(XC*1e-3,YC*1e-3,-ek_depth, 20,cmap=cm.Blues_r)
    plt.colorbar(label="depth (m)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("Ekman Layer depth")
#
    ax2 = plt.subplot(212)
    ax2.set_aspect(1)
    plt.contourf(XC*1e-3,YC*1e-3,mld, 20,cmap=cm.Blues_r)
    plt.colorbar(label="depth (m)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("Mixed Layer depth")
    plt.tight_layout(pad=1)

    if it ==0:
        plt.savefig("./figures/ek_depth_0000"+ str(it) + ".png")
    elif it < 100:
        plt.savefig("./figures/ek_depth_000"+ str(it) + ".png")
    elif it < 1000:
        plt.savefig("./figures/ek_depth_00"+ str(it) + ".png")
    elif it < 10000:
        plt.savefig("./figures/ek_depth_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/ek_depth_"+ str(it) + ".png")
#    plt.close()

"""
    #plot ekman pumping   
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
    if it ==0:
        plt.savefig("./figures/w_mld2_ekmanpump_000"+ str(it) + ".png")
    elif it < 1000:
        plt.savefig("./figures/w_mld2_ekmanpump_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/w_mld2_ekmanpump_"+ str(it) + ".png")
    plt.close()


#plt.contour(xi, yi, zi, v, linewidths=0.5, colors='k')
#plt.contourf(xi, yi, zi, v, cmap=plt.cm.jet)
#x = plt.colorbar(ticks=v)

"""
