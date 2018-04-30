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
Wall=np.ones((40,250,250,120))

plot_abs = 0

#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(240,360,1)
dstart = (itrs[0]+1)*30/720#5
dend = (itrs[-1]+1)*30/720#10
#itrs= [3600]
# time loop
for it in itrs:
    if (it+1)*30==1080:
        W = mit.rdmds('W',((it+2)*30))
    else:
        W = mit.rdmds('W',((it+1)*30))
    Wall[:,:,:,(it-itrs[0])]= W;
    W = []

wmax=np.max(Wall)*1e3
wmin=np.min(Wall)*1e3


#####################################################################

"""
Time averaged
"""
Wmean=np.mean(Wall, axis=3)

wavmax=np.max(Wmean)*1e3
wavmin=np.min(Wmean)*1e3
wavminmax=max(abs(wavmin), abs(wavmax))
v_range = np.linspace(-wavminmax, wavminmax, 101, endpoint=True)
#
#levels = np.linspace(wavmin, wavmax,6, endpoint=True)
levels = np.linspace(-0.02, 0.02,6, endpoint=True)

fig = plt.figure(figsize=(7.5,6))
plt.contourf(XC*1e-3,YC*1e-3,Wmean[25,:,:]*1e3,v_range,cmap=cm.seismic)
#plt.plot(XC[:,166]*1e-3, YC[:,166]*1e-3, ':')
plt.colorbar(label='$\overline{W} \ [mm/s]$', format='%1.3f')

plt.text(10,35,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
plt.text(10,10,'$W_{min}=$ %1.3f $mm/s$' % (wmin))
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title('5 days averaged vertical velocity $\overline{W}$ at $z \sim 230m$, day %d-%d' % (dstart, dend))
#
plt.tight_layout(pad=1)
plt.savefig('./figures/Wmean_200_day%d_%d.png' % (dstart, dend))

#fig = plt.figure(figsize=(14,4))
fig = plt.figure(figsize=(10,6))

#
#ax1 = fig.add_subplot(1, 2, 1)
plt.contourf(YC[:,125]*1e-3,RC.squeeze(),Wmean[:,:,125]*1e3,v_range,cmap=cm.seismic)
#plt.plot(YC[:,125]*1e-3, RC[23].squeeze(), ':')
plt.colorbar(label='$\overline{W} \ [mm/s]$', format='%1.3f')

CS1 = plt.contour(YC[:,125]*1e-3,RC.squeeze(),Wmean[:,:,125]*1e3, levels, colors='k')
plt.clabel(CS1, fmt='%2.3f', colors='k', fontsize=10)

plt.text(10,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
plt.text(10,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title('5 days averaged vertical velocity ($\overline{W}$) at $x=150km$, day %d-%d' % (dstart, dend))

"""
ax2 = fig.add_subplot(1, 2, 2)
plt.contourf(YC[:,166]*1e-3,RC.squeeze(),Wmean[:,:,166]*1e3,v_range,cmap=cm.seismic)
cbar = plt.colorbar(label="W_mean (mm/s)")
#cbar.solids.set_edgecolor("face")
#
CS4 = plt.contour(YC[:,166]*1e-3,RC.squeeze(),Wmean[:,:,166]*1e3, levels)
plt.clabel(CS4, fmt='%2.1f', colors='w', fontsize=10)
#
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title("Time averaged W at x=200km")
"""
plt.tight_layout(pad=1)
plt.savefig('./figures/Wmean_section_day%d_%d.png' % (dstart, dend))
#
#####################################################################
"""
# Mean abs
"""
if plot_abs==1:
    Wmeana=np.mean(abs(Wall), axis=3)
    wavabs=np.max(Wmeana)*1e3
    v_range2 = np.linspace(0, wavabs, 50, endpoint=True)

    fig = plt.figure(figsize=(5,4))
    plt.contourf(XC*1e-3,YC*1e-3,Wmeana[23,:,:]*1e3,v_range2,cmap=cm.Blues)
    plt.plot(XC[:,166]*1e-3, YC[:,166]*1e-3, ':')
    plt.colorbar(label="W_mean (mm/s)")
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("Time averaged W_abs @~200m")
#
    plt.tight_layout(pad=1)
    plt.savefig("./figures/Wmean_abs_200_1.png")

    fig = plt.figure(figsize=(14,4))
#
    ax1 = fig.add_subplot(1, 2, 1)
    plt.contourf(YC[:,125]*1e-3,RC.squeeze(),Wmeana[:,:,125]*1e3,v_range2,cmap=cm.Blues)
    #plt.plot(YC[:,125]*1e-3,RC[23].squeeze(),':')
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
    plt.savefig("./figures/Wmean_abs_xsection_1.png")

"""
"""
###
