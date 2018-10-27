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
timestep = 120
idepth = 50
icx=int(XC[0,:].size/2)
icy=int(YC[0,:].size/2)
plta = 40
pltb = 290

temp_range = np.linspace(0, 1.0, 101, endpoint=True)
co2_range = np.linspace(0, 1.0, 101, endpoint=True)

#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(0,7100,60)
#itrs= [3300,3360, 6000]#, 180, 240, 300,360, 420,480,540,600,660,720,1440,2160,2880,3600]

#fig, ax1 = plt.subplots(figsize=(4,6))
fig = plt.figure(figsize=(12,6))
ax1 = plt.subplot(131)
plt.xlabel(r'$CO_2$ concentration')
plt.ylabel("Depth [m]")
plt.title(r'$CO_2$ concentration at initial eddy core')
plt.grid(True)
ax2 = plt.subplot(132,  sharex=ax1, sharey=ax1)
plt.xlabel(r'$CO_2$ concentration')
#plt.ylabel("Depth [m]")
plt.title(r'$CO_2$ at %d km north from initial eddy core' %((XC[100,100]-XC[90,90])*1e-3))
plt.grid(True)
ax3 = plt.subplot(133,  sharex=ax1, sharey=ax1)
plt.xlabel(r'$CO_2$ concentration')
#plt.ylabel("Depth [m]")
plt.title(r'$CO_2$ at %d km south from initial eddy core' %((XC[80,80]-XC[90,90])*1e-3))
plt.grid(True)

# time loop
for it in itrs:
    co2 = mit.rdmds('PTRACER01',it)
#    visc = mit.rdmds('KPPviscAz',it)
    
# Surface
    plt.figure(figsize=(5,4))
    plt.contourf(XC[plta:pltb,plta:pltb]*1e-3,YC[plta:pltb,plta:pltb]*1e-3,co2[0,plta:pltb,plta:pltb],temp_range, cmap=cm.terrain_r)
    plt.colorbar(label='concentration')
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.text(400,400,r'$CO_2$ max = %.2f' % (np.max(co2)))
    plt.title(r'Surface $CO_2$ @ timestep = ' + str(int(it*timestep/3600))+ 'hr')
    plt.tight_layout(pad=1)
    if it==0:
        plt.savefig("./figures/surf_co2_0000"+ str(it) + ".png")
    elif it<100:
        plt.savefig("./figures/surf_co2_000"+ str(it) + ".png")
    elif it<1000:
        plt.savefig("./figures/surf_co2_00"+ str(it) + ".png")
    elif it<10000:
        plt.savefig("./figures/surf_co2_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/surf_co2_"+ str(it) + ".png")
    plt.close()

# Section
    fig = plt.figure(figsize=(7,4))
#
#    ax5 = fig.add_subplot(1, 2, 1)
    plt.contourf(YC[plta:pltb, icx]*1e-3,RC[:idepth].squeeze(),co2[:idepth,plta:pltb,icx],co2_range, cmap=cm.terrain_r)
    plt.colorbar(label="concentration")
    plt.xlabel("y (km)")
    plt.ylabel("z (m)")
    plt.title(r'$CO_2$ @ x=150km, timestep = ' + str(int(it*timestep/3600)) + 'hr')
    """
    ax6 = fig.add_subplot(1, 2, 2)
    plt.contourf(XC[icy,plta:pltb]*1e-3,RC[:idepth].squeeze(),co2[:idepth,icy,plta:pltb],co2_range, cmap=cm.rainbow)
    plt.colorbar(label="concentration")
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title(r'$CO_2$ @ y=150km, timestep = ' + str(it))
    """
    plt.tight_layout(pad=1)
    if it==0:
        plt.savefig("./figures/co2_xsection_0000"+ str(it) + ".png", dpi=300)
    elif it<100:
        plt.savefig("./figures/co2_xsection_000"+ str(it) + ".png", dpi=300)
    elif it<1000:
        plt.savefig("./figures/co2_xsection_00"+ str(it) + ".png", dpi=300)
    elif it<10000:
        plt.savefig("./figures/co2_xsection_0"+ str(it) + ".png", dpi=300)
    else:
        plt.savefig("./figures/co2_xsection_"+ str(it) + ".png", dpi=300)
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
    if (it % 720) == 0:
#    fig = plt.figure(figsize=(5,7))
        ax1.plot(co2[:,icy,icx],RC.squeeze())
        ax2.plot(co2[:,100,90],RC.squeeze(), label='time step %d day' % (it*timestep/86400))
        ax3.plot(co2[:,80,90],RC.squeeze())

ax2.legend(loc='lower right')
plt.tight_layout(pad=1)
plt.savefig('./figures/co2_profile.png')

plt.ylim(-100,0)
plt.tight_layout(pad=1)
plt.savefig('./figures/co2_profile100.png')

plt.ylim(-300,0)
plt.tight_layout(pad=1)
plt.savefig('./figures/co2_profile300.png')

