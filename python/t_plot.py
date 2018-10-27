"""
Plot the Temperature on the surface and the cross section
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm

plt.ion()

# Set default fontsizes for plots
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')
timestep = 180
idepth = 65
icx=int(XC[0,:].size/2)
icy=int(YC[0,:].size/2)

temp_range = np.linspace(20, 25, 101, endpoint=True)

#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
#itrs=np.arange(1170,10980,180)
itrs= [0, 720, 3600,7200,10800, 14400]

#fig, ax1 = plt.subplots(figsize=(4,6))
fig = plt.figure(figsize=(12,6))
ax1 = plt.subplot(131)
plt.xlabel('Temperature')
plt.ylabel("Depth [m]")
plt.title('Temperature at initial eddy core')
plt.grid(True)
ax2 = plt.subplot(132,  sharex=ax1, sharey=ax1)
plt.xlabel('Temperature')
#plt.ylabel("Depth [m]")
plt.title('T at %d km north from initial eddy core' %((XC[100,100]-XC[90,90])*1e-3))
plt.grid(True)
ax3 = plt.subplot(133,  sharex=ax1, sharey=ax1)
plt.xlabel('Temperature')
#plt.ylabel("Depth [m]")
plt.title('T at %d km south from initial eddy core' %((XC[80,80]-XC[90,90])*1e-3))
plt.grid(True)

# time loop
for it in itrs:
    Temp = mit.rdmds('T',it)
#    visc = mit.rdmds('KPPviscAz',it)
    
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
        plt.savefig("./figures/SST_0000"+ str(it) + ".png")
    elif it<1000:
        plt.savefig("./figures/SST_00"+ str(it) + ".png")
    elif it<10000:
        plt.savefig("./figures/SST_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/SST_"+ str(it) + ".png")
#    plt.close()

# Section
    fig = plt.figure(figsize=(7,4))
#
#    ax1 = fig.add_subplot(1, 2, 1)
    plt.contourf(XC[icy,:]*1e-3,RC.squeeze(),Temp[:,:,icx],100, cmap=cm.nipy_spectral)
    plt.colorbar(label="Temp (°C)")
    plt.xlabel("y (km)")
    plt.ylabel("z (m)")
    plt.title("Temp @ x=150km, timestep = " + str(it))
    """
    ax2 = fig.add_subplot(1, 2, 2)
    plt.contourf(YC[:,icx]*1e-3,RC.squeeze(),Temp[:,icy,:],100, cmap=cm.nipy_spectral)
    plt.colorbar(label="Temp (°C)")
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("Temp @ y=150km, timestep = " + str(it))
    """
    plt.tight_layout(pad=1)
    if it==0:
        plt.savefig("./figures/T_section_0000"+ str(it) + ".png", dpi=300)
    elif it<1000:
        plt.savefig("./figures/T_section_00"+ str(it) + ".png", dpi=300)
    else:
        plt.savefig("./figures/T_section_"+ str(it) + ".png", dpi=300)
    plt.close()
    """
    if it==0:
        plt.savefig("./figures/T_section_0000"+ str(it) + ".eps", format='eps', dpi=300)
    elif it<1000:
        plt.savefig("./figures/T_section_00"+ str(it) + ".eps", format='eps', dpi=300)
    else:
        plt.savefig("./figures/T_section_"+ str(it) + ".eps", format='eps', dpi=300)
#    plt.close()
    """
#
#    fig = plt.figure(figsize=(5,7))
    ax1.plot(Temp[0:idepth,icy,icx],RC[0:idepth].squeeze())
    ax2.plot(Temp[0:idepth,100,90],RC[0:idepth].squeeze(), label='time step %d day' % (it*timestep/86400))
    ax3.plot(Temp[0:idepth,80,90],RC[0:idepth].squeeze())

ax2.legend(loc='lower right')
plt.tight_layout(pad=1)
plt.savefig('./figures/T_profile.png')

plt.ylim(-100,0)
plt.tight_layout(pad=1)
plt.savefig('./figures/T_profile100.png')

plt.ylim(-300,0)
plt.tight_layout(pad=1)
plt.savefig('./figures/T_profile300.png')

