"""
Plot the w on the surface and the cross section
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

plt.ion()

timestep = 120 #timestep input in second
dumpfreq = 7200 # file every 3 hours for dumpfreq 10800
file_dit = dumpfreq/timestep
day_s = 0
day_e = 15
startfile = day_s*file_dit #or any integer*file_dit
endfile = day_e*86400/timestep + 50 #1 day 
itrs = np.arange(startfile,endfile, file_dit)
itrs= [60, 120] #2400, 3600, 4800,7200]#
itrs = mit.mds.scanforfiles('diagTAUX')

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
#RC = mit.rdmds('RC')

#v_range = 100
v_range_max = 0.2
v_range = np.linspace(0.1, v_range_max, 101, endpoint=True) #.02 or .12 
v_ticks = np.linspace(0.1, v_range_max, 11, endpoint=True)
#levels = np.linspace(-0.02, 0.02,6, endpoint=True)
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)

Rmax = 25*1e3 #25 15.625*1e3
pos_tx = 400
pos_ty = 400 #440 500
nx = XC.shape[1]
ny = YC.shape[0]

icx = int(nx/2) 
icy = int(ny/2)

plta = 40
pltb = 290

xc_dom =XC[plta:pltb, plta:pltb]*1e-3
yc_dom =YC[plta:pltb, plta:pltb]*1e-3

nu=0.01 #from kpp=0.01 or Az=4*1e-3
f0=6.5e-5
de = -np.sqrt(2*nu/f0)

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

# time loop
for it in itrs:
    taux = mit.rdmds('diagTAUX',it)
    tauy = mit.rdmds('diagTAUY',it)
#
    """
    """
# ===========================================================================================
# plot vertically averaged W and section

    fig = plt.figure(figsize=(5.5,9))
    ax1 = plt.subplot(211)
    ax1.set_aspect(1)
    plt.contourf(xc_dom,yc_dom,taux[plta:pltb,plta:pltb],v_range,cmap=cm.rainbow)#
    plt.plot(XC[plta:pltb,icx]*1e-3,YC[plta:pltb,icx]*1e-3, ':')#
    cb=plt.colorbar(format='%1.2f', ticks=v_ticks)
    cb.set_label(r'$\tau \, (N/m^2)$', labelpad=-40, y=1.1, rotation=0)
    cb.ax.tick_params(labelsize=10)
    CS2 = plt.contour(xc_dom,yc_dom,taux[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
    plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
    rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
    ax1.add_artist(rmax1)
    rmax2 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
    ax1.add_artist(rmax2)
    plt.text(pos_tx,pos_ty,r'$\tau_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (np.min(taux), np.max(taux)), fontsize=12)
#
    plt.xlabel("x (km)", fontsize=10)
    plt.ylabel("y (km)", fontsize=10)
#
    plt.title(r'Windstress, $\tau_x$, timestep %d hr' % (int(it*timestep/3600)),fontsize=11)

#
    ax2 = plt.subplot(212, sharex=ax1)
    plt.contourf(xc_dom,yc_dom,tauy[plta:pltb,plta:pltb],v_range,cmap=cm.rainbow)#
    plt.plot(XC[plta:pltb,icx]*1e-3,YC[plta:pltb,icx]*1e-3, ':')#
    cb=plt.colorbar(format='%1.2f', ticks=v_ticks)
    cb.set_label(r'$\tau \, (N/m^2)$', labelpad=-40, y=1.1, rotation=0)
    cb.ax.tick_params(labelsize=10)
    CS2 = plt.contour(xc_dom,yc_dom,tauy[plta:pltb,plta:pltb], levels, colors='0.6')
    plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
    rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
    ax1.add_artist(rmax1)
    rmax2 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
    ax1.add_artist(rmax2)
    plt.text(pos_tx,pos_ty,r'$\tau_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (np.min(tauy), np.max(tauy)), fontsize=12)
#
    plt.xlabel("x (km)", fontsize=10)
    plt.ylabel("y (km)", fontsize=10)
#
    plt.title(r'Windstress, $\tau_y$, timestep %d hr' % (int(it*timestep/3600)),fontsize=11)
#
    if it == 0:
        plt.savefig('./figures/tau_0000%d.png' % (it))
    elif it<100:
        plt.savefig('./figures/tau_000%d.png' % (it))
    elif it<1000:
        plt.savefig('./figures/tau_00%d.png' % (it))
    elif it<10000:
        plt.savefig('./figures/tau_0%d.png' % (it))
    else:
        plt.savefig('./figures/tau_%d.png' % (it))
#    plt.close()


