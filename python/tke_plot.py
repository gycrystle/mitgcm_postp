"""
Plot the TKE on the surface and the cross section
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

v_range = np.linspace(-0.6, 0.6, 101, endpoint=True)
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
#itrs=np.arange(1170,10980,180)
itrs= [2160]

# time loop
for it in itrs:
    U = mit.rdmds('U',it)
    V = mit.rdmds('V',it)
    TKE = mit.rdmds('GGL90TKE', it)
    umax = np.max(U)
    vmax = np.max(V)
    #fig = plt.figure(figsize=(10,4))
    fig = plt.figure(figsize=(15,6))
#
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    plt.contourf(XC*1e-3,YC*1e-3,TKE[1,:,:],100,cmap=cm.seismic)
    #plt.plot(XC[:,166]*1e-3,YC[:,166]*1e-3, ':')
    plt.colorbar(label='$U \ [m/s]$', format='%1.3f')
#    CS1 = plt.contour(XC[125,:]*1e-3,YC[:,125]*1e-3,U[1,:,:], levels, colors='k')
#    plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)

    #plt.text(10,-1250,'$U_{max}=$ %1.3f $m/s$' % (wmax))
#    plt.text(10,-1350,'$U_{max}=$ %1.3f $m/s$' % (umax))

    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title('$TKE_{surf}$, timestep = ' + str(it))
#
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_aspect(1)
    plt.contourf(XC*1e-3,YC*1e-3,TKE[1,:,:],v_range,cmap=cm.seismic)#
    plt.colorbar(label='$V \ [m/s]$', format='%1.3f')
    
#    CS2 = plt.contour(XC[125,:]*1e-3,YC[:,125]*1e-3,V[1,:,:], levels, colors='k')
#    plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)

    #plt.text(10,-1250,'$U_{max}=$ %1.3f $m/s$' % (wmax))
#    plt.text(10,-1350,'$V_{max}=$ %1.3f $m/s$' % (vmax))

    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title('$V_{surf}$, timestep = ' + str(it))

    plt.tight_layout(pad=1)
    if it==0:
        plt.savefig("./figures/TKE_surf_0000"+ str(it) + ".png")
    elif it<100:
        plt.savefig("./figures/TKE_surf_000"+ str(it) + ".png")
    elif it<1000:
        plt.savefig("./figures/TKE_surf_00"+ str(it) + ".png")
    elif it<10000:
        plt.savefig("./figures/TKE_surf_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/TKE_surf_"+ str(it) + ".png")
#    plt.close()

    fig = plt.figure(figsize=(10,6))
#
    #ax1 = fig.add_subplot(1, 2, 1)
    plt.contourf(YC[:,125]*1e-3,RC.squeeze(),TKE[:,:,125],100,cmap=cm.seismic)
    plt.colorbar(label='$U \ [m/s]$', format='%1.3f')
#    
#    CS3 = plt.contour(YC[:,125]*1e-3,RC.squeeze(),U[:,:,125], levels, colors='k')
#    plt.clabel(CS3, fmt='%2.2f', colors='k', fontsize=10)

#    plt.text(10,-1350,'$U_{max}=$ %1.3f $m/s$' % (umax))

    plt.xlabel("y (km)")
    plt.ylabel("z (m)")
    plt.title("$U$ at $x=150km$, timestep = " + str(it))
#
    """
    ax2 = fig.add_subplot(1, 2, 2)
    plt.contourf(XC[125,:]*1e-3,RC.squeeze(),V[:,125,:],v_range,cmap=cm.seismic)
    #plt.contourf(YC[:,166]*1e-3,RC.squeeze(),U[:,:,166],v_range,cmap=cm.seismic)
    plt.colorbar(label='$V \ [m/s]$', format='%1.3f')

    CS4 = plt.contour(XC[125,:]*1e-3,RC.squeeze(),V[:,125,:], levels, colors='k')
    plt.clabel(CS4, fmt='%2.2f', colors='k', fontsize=10)

    #plt.text(10,-1250,'$U_{max}=$ %1.3f $m/s$' % (wmax))
    plt.text(10,-1350,'$V_{max}=$ %1.3f $m/s$' % (vmax))

    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title('$V$ @ $y=150km$, timestep = ' + str(it))
    """
    plt.tight_layout(pad=1)
    if it==0:
        plt.savefig("./figures/TKE_section_0000"+ str(it) + ".png")
    elif it<100:
        plt.savefig("./figures/TKE_section_000"+ str(it) + ".png")
    elif it<1000:
        plt.savefig("./figures/TKE_section_00"+ str(it) + ".png")
    elif it<10000:
        plt.savefig("./figures/TKE_section_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/TKE_section_"+ str(it) + ".png")
#    plt.close()
###
