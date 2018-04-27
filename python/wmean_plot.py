"""
Plot the time average of w
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
Wall=np.ones((40,250,250,60))

v_range = np.linspace(-0.1, 0.1, 101, endpoint=True)
v_range2 = np.linspace(0, 0.2, 51, endpoint=True)


#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(0,60,1)
#itrs= [3600]
# time loop
for it in itrs:
    W = mit.rdmds('W',((it+1)*60))
    Wall[:,:,:,it]=W
#    W=[]
#    V = mit.rdmds('V',it)
"""
Time averaged
"""
Wmean=np.mean(Wall, axis=3)
fig = plt.figure(figsize=(5,4))
plt.contourf(XC*1e-3,YC*1e-3,Wmean[12,:,:]*1e3,v_range,cmap=cm.seismic)
plt.plot(XC[:,166]*1e-3, YC[:,166]*1e-3, ':')
plt.colorbar(label="W_mean (mm/s)")
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Time averaged W @125m")
#
plt.tight_layout(pad=1)
plt.savefig("./figures/Wmean_125.png")

fig = plt.figure(figsize=(14,4))
#
ax1 = fig.add_subplot(1, 2, 1)
plt.contourf(YC[:,125]*1e-3,RC.squeeze(),Wmean[:,:,125]*1e3,v_range,cmap=cm.seismic)
plt.colorbar(label="W_mean (m/s)")
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Time averaged W at x=150km")

ax2 = fig.add_subplot(1, 2, 2)
plt.contourf(YC[:,166]*1e-3,RC.squeeze(),Wmean[:,:,166]*1e3,v_range,cmap=cm.seismic)
plt.colorbar(label="W_mean (m/s)")
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Time averaged W at x=200km")

plt.tight_layout(pad=1)
plt.savefig("./figures/Wmean_section.png")
#

"""
# Mean abs
"""
Wmeana=np.mean(abs(Wall), axis=3)
fig = plt.figure(figsize=(5,4))
plt.contourf(XC*1e-3,YC*1e-3,Wmeana[12,:,:]*1e3,100,cmap=cm.Blues)
plt.plot(XC[:,166]*1e-3, YC[:,166]*1e-3, ':')
plt.colorbar(label="W_mean (mm/s)")
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Time averaged W_abs @125m")
#
plt.tight_layout(pad=1)
plt.savefig("./figures/Wmean_abs_125.png")

fig = plt.figure(figsize=(14,4))
#
ax1 = fig.add_subplot(1, 2, 1)
plt.contourf(YC[:,125]*1e-3,RC.squeeze(),Wmeana[:,:,125]*1e3,v_range2,cmap=cm.Blues)
#plt.plot(YC[:,125]*1e-3,RC[22].squeeze(),':')
plt.colorbar(label="W_mean (m/s)")
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Time averaged W_abs @ x=150km")

ax2 = fig.add_subplot(1, 2, 2)
plt.contourf(YC[:,166]*1e-3,RC.squeeze(),Wmeana[:,:,166]*1e3,v_range2,cmap=cm.Blues)
plt.colorbar(label="W_mean (m/s)")
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Time averaged W_abs @ x=200km")

plt.tight_layout(pad=1)
plt.savefig("./figures/Wmean_abs_xsectionio.png")

"""
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_aspect(1)
    plt.contourf(XC*1e-3,YC*1e-3,V[1,:,:],v_range,cmap=cm.seismic)#
    plt.colorbar(label="V (m/s)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("V_surf @ timestep = " + str(it))

    plt.tight_layout(pad=1)
    plt.savefig("./figures/V_surf_"+ str(it) + ".png")

    fig = plt.figure(figsize=(10,4))
#
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    plt.contourf(XC*1e-3,YC*1e-3,U[1,:,:],v_range,cmap=cm.seismic)
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
    plt.savefig("./figures/V_surf_"+ str(it) + ".png")

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
    plt.contourf(XC[125,:]*1e-3,RC.squeeze(),V[:,125,:],v_range,cmap=cm.seismic)
    plt.colorbar(label="V (m/s)")
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("V @ x=150km, timestep = " + str(it))

    plt.tight_layout(pad=1)
    plt.savefig("./figures/V_section_"+ str(it) + ".png")
"""
###
