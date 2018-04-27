"""
Plot the Temperature on the surface and the cross section
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

temp_range = np.linspace(26, 30, 101, endpoint=True)

#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(0,3780,180)
#itrs= [3600]

# time loop
for it in itrs:
    Temp = mit.rdmds('T',it)
    
# Temp plot
# SST
    plt.figure(figsize=(5,4))
    plt.contourf(XC*1e-3,YC*1e-3,Temp[1,:,:],temp_range, cmap=cm.rainbow)
    plt.colorbar(label="Temp (°C)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("SST @ timestep = " + str(it))
    plt.tight_layout(pad=1)
    if it==0:
        plt.savefig("./figures/SST_000"+ str(it) + ".png")
    elif it<1000:
        plt.savefig("./figures/SST_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/SST_"+ str(it) + ".png")
    plt.close()

# Section
    fig = plt.figure(figsize=(14,4))
#
    ax1 = fig.add_subplot(1, 2, 1)
    plt.contourf(XC[125,:]*1e-3,RC.squeeze(),Temp[:,:,125],100, cmap=cm.nipy_spectral)
    plt.colorbar(label="Temp (°C)")
    plt.xlabel("y (km)")
    plt.ylabel("z (m)")
    plt.title("Temp @ x=150km, timestep = " + str(it))

    ax2 = fig.add_subplot(1, 2, 2)
    plt.contourf(YC[:,125]*1e-3,RC.squeeze(),Temp[:,125,:],100, cmap=cm.nipy_spectral)
    plt.colorbar(label="Temp (°C)")
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("Temp @ y=150km, timestep = " + str(it))

    plt.tight_layout(pad=1)
    if it==0:
        plt.savefig("./figures/T_section_000"+ str(it) + ".png")
    elif it<1000:
        plt.savefig("./figures/T_section_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/T_section_"+ str(it) + ".png")
    plt.close()


