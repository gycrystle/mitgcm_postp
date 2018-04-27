"""
Plot the u, v on the surface and the cross section
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm

#plt.ion()

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')

v_range = np.linspace(-0.6, 0.6, 101, endpoint=True)


#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(0,3780,180)
#itrs= [3600]

# time loop
for it in itrs:
    U = mit.rdmds('U',it)
    V = mit.rdmds('V',it)

    fig = plt.figure(figsize=(10,4))
#
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    plt.contourf(XC*1e-3,YC*1e-3,U[1,:,:],v_range,cmap=cm.seismic)
    plt.plot(XC[:,166]*1e-3,YC[:,166]*1e-3, ':')
    plt.colorbar(label="U (m/s)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("U_surf @ timestep = " + str(it))
#
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_aspect(1)
    plt.contourf(XC*1e-3,YC*1e-3,V[1,:,:],v_range,cmap=cm.seismic)#
    plt.colorbar(label="V (m/s)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("V_surf @ timestep = " + str(it))

    plt.tight_layout(pad=1)
    if it==0:
        plt.savefig("./figures/U_surf_000"+ str(it) + ".png")
    elif it<100:
        plt.savefig("./figures/U_surf_00"+ str(it) + ".png")
    elif it<1000:
        plt.savefig("./figures/U_surf_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/U_surf_"+ str(it) + ".png")
    plt.close()

    fig = plt.figure(figsize=(14,4))
#
    ax1 = fig.add_subplot(1, 2, 1)
    plt.contourf(YC[:,125]*1e-3,RC.squeeze(),U[:,:,125],v_range,cmap=cm.seismic)
    plt.colorbar(label="U (m/s)")
    plt.xlabel("y (km)")
    plt.ylabel("z (m)")
    plt.title("U @ x=150km, timestep = " + str(it))
#
    ax2 = fig.add_subplot(1, 2, 2)
    plt.contourf(YC[:,166]*1e-3,RC.squeeze(),U[:,:,166],v_range,cmap=cm.seismic)
    plt.colorbar(label="U (m/s)")
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("U @ x=200km, timestep = " + str(it))

    plt.tight_layout(pad=1)
    if it==0:
        plt.savefig("./figures/U_section_000"+ str(it) + ".png")
    elif it<100:
        plt.savefig("./figures/U_section_00"+ str(it) + ".png")
    elif it<1000:
        plt.savefig("./figures/U_section_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/U_section_"+ str(it) + ".png")
    plt.close()
###
