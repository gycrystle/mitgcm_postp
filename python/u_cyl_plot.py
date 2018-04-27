"""
Plot the u, v on the surface and the cross section
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
theta = np.arctan2(YC-150000, XC-150000)
#R = np.hypot(XC,YC)
v_range = np.linspace(-0.5, 0.5, 101, endpoint=True)


#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
#itrs=np.arange(0,3780,180)
itrs= [2160, 3600]

# time loop
for it in itrs:
    U = mit.rdmds('U',it)
    V = mit.rdmds('V',it)
    U_theta = (V * np.cos(theta)) - (U * np.sin(theta));


    fig = plt.figure(figsize=(5,4))
#
#    ax1 = fig.add_subplot(1, 2, 1)
#    ax1.set_aspect(1)
    plt.contourf(XC*1e-3,YC*1e-3,U_theta[1,:,:],v_range,cmap=cm.seismic)
    plt.colorbar(label="U_theta (m/s)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("U_theta surface @ timestep = " + str(it))
#

    plt.tight_layout(pad=1)
    plt.savefig("./figures/U_theta_surf_"+ str(it) + ".png")


    fig = plt.figure(figsize=(7,4))
#
#    ax1 = fig.add_subplot(1, 2, 1)
#    plt.contourf(YC[:,125]*1e-3,RC.squeeze(),U[:,:,125],v_range,cmap=cm.seismic)
#    plt.colorbar(label="U (m/s)")
#    plt.xlabel("y (km)")
#    plt.ylabel("z (m)")
#    plt.title("U @ x=150km, timestep = " + str(it))
#
#    ax2 = fig.add_subplot(1, 2, 2)
    plt.contourf(XC[125,:]*1e-3,RC.squeeze(),U_theta[:,125,:],v_range,cmap=cm.seismic)
    plt.colorbar(label="U_theta (m/s)")
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("U_theta @ x=150km, timestep = " + str(it))

    plt.tight_layout(pad=1)
    plt.savefig("./figures/U_theta_section_"+ str(it) + ".png")

    
###
