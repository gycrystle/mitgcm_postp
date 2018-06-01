"""
Python script to read compute and plot Vorticity (zeta) output from MITgcm
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)
    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))


plt.ion()

nit = 60 # no of time step 
ts = 60  # time step = ts*dt (in second); 
nr = 100  # no of grid vertical 
nx = 180 # no of grid in x
ny = 180 # no of grid in y
sec_y = int(nx/2)

#vort_range = np.linspace(-0.5, 0.5, 101, endpoint=True)
vort_range = 101

itrs = np.arange(0, 7210, 120)
#itrs = [2160, 2880, 3600] #, 2160, 2880, 3600, 4320, 5040, 5760, 6480, 7200]


rho0 = 999.8
XC = mit.rdmds('XC')
dXC = mit.rdmds('DXC')
YC = mit.rdmds('YC')
dYC = mit.rdmds('DYC')
RC = mit.rdmds('RC')

dyU = np.zeros((nr,ny,nx));
dxV = np.zeros((nr,ny,nx));
         
zeta = np.zeros((nr,ny,nx));

for it in itrs:
    U = mit.rdmds('U',it)
    V = mit.rdmds('V',it)
#
    for i in range(1,nx):
        for j in range(1,ny):
            #for k in range(0,nr):
            dyU[:,j,i] = (U[:,j,i]-U[:,j-1,i])/dYC[j,i];
            dxV[:,j,i] = (V[:,j,i]-V[:,j,i-1])/dXC[j,i];
            zeta[:,j,i] = dxV[:,j,i]-dyU[:,j,i];
#
    zetamin = np.min(zeta)*1e5
#############Plot figures ##############
    #plt.figure()
    fig = plt.figure(figsize=(12,4))
#
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    plt.contourf(XC*1e-3,YC*1e-3,zeta[0,:,:]*1e5, 100,vmin=zetamin,vmax=-zetamin,cmap=cm.seismic)
    plt.colorbar(label=r'$\zeta \ 10^5 s^{-1}$')
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(r'Surface vorticity, $\zeta \ 10^5 s^{-1}$')
#
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_aspect(2)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf (XC[sec_y,:]*1e-3, RC.squeeze(), zeta[:, sec_y,:]*1e5, 100,vmin=zetamin,vmax=-zetamin,cmap=cm.seismic)
    plt.colorbar(label=r'$\zeta \ 10^5 s^{-1}$')
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("Vorticity")

    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/vort_000"+ str(it) + ".png")
    elif it < 1000:
        plt.savefig("./figures/vort_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/vort_"+ str(it) + ".png")
    plt.close()
"""
##========== from diag ====================================
itrsd = np.arange(0, 3600, 720)

for it in itrsd:
    vort = mit.rdmds('momvort',it)
    KE = mit.rdmds('momKE',it)

    fig = plt.figure(figsize=(12,4))    

    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf(XC*1e-3,YC*1e-3,vort[0,:,:]*1e5, vort_range,MidpointNormalize(vmax=-vmin, midpoint=0.0),cmap=cm.seismic)
    plt.colorbar(label=r'$\zeta \ 10^5 s^{-1}$')
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(r'Surface vorticity, $\zeta \ 10^5 s^{-1}$')
    #plt.savefig("mld"+ str(i) + ".png")

    ax7 = fig.add_subplot(1, 2, 2)
    ax7.set_aspect(0.01)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf(XC[sec_y,:]*1e-3,RC.squeeze(),vort[:,sec_y,:]*1e5, vort_range,cmap=cm.seismic)
    plt.colorbar(label=r'$\zeta \ 10^5 s^{-1}$')
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("Vorticity")

    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/vortmom_000"+ str(it) + ".png")
    elif it < 1000:
        plt.savefig("./figures/vortmom_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/vortmom_"+ str(it) + ".png")
#    plt.close()

# KE
    fig = plt.figure(figsize=(12,4))    

    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf(XC*1e-3,YC*1e-3,KE[0,:,:]*1e5, 100,cmap=cm.inferno)
    plt.colorbar(label='KE')
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("KE")
    #plt.savefig("mld"+ str(i) + ".png")

    ax7 = fig.add_subplot(1, 2, 2)
    ax7.set_aspect(0.1)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf(XC[sec_y,:]*1e-3,RC.squeeze(),KE[:,sec_y,:]*1e5, 100,cmap=cm.hsv)
    plt.colorbar(label=r'$\zeta$')
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("KE")

    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/KE_000"+ str(it) + ".png")
    elif it < 1000:
        plt.savefig("./figures/KE_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/KE_"+ str(it) + ".png")
#    plt.close()
"""

#plt.contour(xi, yi, zi, v, linewidths=0.5, colors='k')
#plt.contourf(xi, yi, zi, v, cmap=plt.cm.jet)
#x = plt.colorbar(ticks=v)

