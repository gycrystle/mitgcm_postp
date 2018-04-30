"""
Plot the w on the surface and the cross section
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm

plt.ion()

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')

v_range = np.linspace(-0.4, 0.4, 101, endpoint=True)
#levels = np.linspace(-0.02, 0.02,6, endpoint=True)
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(1170,10890,90)
#itrs= [2160, 3600]

# time loop
for it in itrs:
    W = mit.rdmds('W',it)
    wmax=np.max(W)*1e3
    wmin=np.min(W)*1e3

#
    fig = plt.figure(figsize=(7.5,6))
#
#    ax1 = fig.add_subplot(1, 2, 1)
#    ax1.set_aspect(1)
    plt.contourf(XC*1e-3,YC*1e-3,W[25,:,:]*1e3,v_range,cmap=cm.seismic)#
    #plt.plot(XC[:,166]*1e-3,YC[:,166]*1e-3, ':')#axvline(x=XC[166,166])
    plt.colorbar(label="W (mm/s)")
    CS2 = plt.contour(XC[125,:]*1e-3,YC[:,125]*1e-3,W[25,:,:]*1e3, levels, colors='0.6')
    plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)

    plt.text(10,35,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
    plt.text(10,10,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title('Vertical velocity $W$ at $z \sim 230m$, timestep %d' % (it))
    if it==0:
        plt.savefig("./figures/W_230_000"+ str(it) + ".png")
    elif it<100:
        plt.savefig("./figures/W_230_00"+ str(it) + ".png")
    elif it<1000:
        plt.savefig("./figures/W_230_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/W_230_"+ str(it) + ".png")
    plt.close()
#
    fig = plt.figure(figsize=(10,6))
#
    #ax1 = fig.add_subplot(1, 2, 1)
    plt.contourf(YC[:,125]*1e-3,RC.squeeze(),W[:,:,125]*1e3,v_range,cmap=cm.seismic)
    plt.colorbar(label="W (mm/s)")
    CS1 = plt.contour(YC[:,125]*1e-3,RC.squeeze(),W[:,:,125]*1e3, levels, colors='0.6')
    plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)

    plt.text(10,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
    plt.text(10,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("y (km)")
    plt.ylabel("z (m)")
    plt.title('Vertical velocity $W$ at $x=150km$, timestep %d' % (it))
#    plt.savefig("./figures/U_surf"+ str(it) + ".png")

    #
    """
    ax2 = fig.add_subplot(1, 2, 2)
    plt.contourf(YC[:,166]*1e-3,RC.squeeze(),W[:,:,166]*1e3,v_range,cmap=cm.seismic)
    plt.colorbar(label="W (mm/s)")
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("W @ x=200km, timestep = " + str(it))
    """
    plt.tight_layout(pad=1)
    if it==0:
        plt.savefig("./figures/W_section_000"+ str(it) + ".png")
    elif it<100:
        plt.savefig("./figures/W_section_00"+ str(it) + ".png")
    elif it<1000:
        plt.savefig("./figures/W_section_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/W_section_"+ str(it) + ".png")
    plt.close()

