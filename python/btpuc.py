"""
Plot the w_ekman from ek transport divergence
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors

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


class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)
    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

# select plot domain
plta = 30
pltb = 150
idepth = 16
depth_plot_ratio = 0.7

timestep = 180 #timestep input in second
dumpfreq = 10800 # in second

day_s = 0
day_e = 5

ts = dumpfreq/timestep  # time step = ts*dt (in second); = 7200 = dumpfreq

nstart = int(day_s*86400/dumpfreq) # integer of dumpfrec, either start from zero or 
nend = int(day_e*86400/dumpfreq) # no of time step 
itpd = int(86400/dumpfreq)

f0 = 1e-4
rho0 = 999.8
Rmax = 25e3
vmax = 0.4

#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(nstart,nend+1,1)
#itrs=[0,1,20,120]#480,960,1440,1920,2400,7200]#,480,960,1440,1920,2400]
nit=itrs.size
#nit=4

vr_range=np.linspace(-0.1, 0.1, 100, endpoint=True)
vr_ticks = np.linspace(-0.1, 0.1, 11, endpoint=True)
vt_range=np.linspace(-0.45, 0, 100, endpoint=True)
vt_ticks = np.linspace(-0.45, 0, 11, endpoint=True)

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')
XG = mit.rdmds('XG')
YG = mit.rdmds('YG')

# indices of eddy center
icx = int(XC[0,:].size/2)
icy = int(YC[:,0].size/2)

#transform to polar coord
thetaC = np.arctan2(YC-YC[icy, icx], XC-XC[icy, icx])
radC = np.hypot(XC-XC[icy, icx],YC-YC[icy, icx])

nr = RC.size
nx = XC[0,:].size
ny = YC[:,0].size

U_c = np.zeros((nr,ny,nx))
V_c = np.zeros((nr,ny,nx))

xc_dom =XC[plta:pltb, plta:pltb]*1e-3
yc_dom =YC[plta:pltb, plta:pltb]*1e-3

dstart = day_s #5
dend = day_e # (itrs[-1]+1)*ts/720  #10

levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
# time loop start
for it in itrs:
    U = mit.rdmds('U',(it*ts))
    V = mit.rdmds('V',(it*ts))
# Must interpolate to cell center before computing the tangential and radial velocity
    for z in np.arange(0,nr):
        for y in np.arange(0,ny):
            U_c[z,y,:] = np.interp(XC[y,:],XG[y,:],U[z,y,:])
        for x in np.arange(0,nx):
            V_c[z,:,x] = np.interp(YC[:,x],YG[:,x],V[z,:,x])

    U_theta = V_c*np.cos(thetaC) - U_c*np.sin(thetaC)
    U_r = U_c*np.cos(thetaC)+V_c*np.sin(thetaC)
#
    it_hour = int(it*dumpfreq/3600)
    fig = plt.figure(figsize=(5.5,9))
#
    ax1 = plt.subplot(211)
    ax1.set_aspect(1)
    plt.contourf(xc_dom,yc_dom,U_theta[idepth,plta:pltb,plta:pltb],vt_range,cmap=cm.Blues_r)#
    cb=plt.colorbar(format='%1.3f', ticks=vt_ticks)
    cb.set_label(r'$U_\theta \, (m/s)$', labelpad=-40, y=1.1, rotation=0)
    cl = plt.contour(xc_dom,yc_dom,U_theta[idepth,plta:pltb,plta:pltb], levels, colors='0.6')
    plt.clabel(cl, fmt='%2.2f', colors='k', fontsize=10)
    rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
    ax1.add_artist(rmax1)
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(r'Tangential velocity, $U_\theta \, (z=%d m)$, timestep %d hr' % (RC[idepth], it_hour),fontsize=11)
#plt.text(480,480,'$U_{extrema}=$ [%1.3f, %1.3f] $m/s$' % (np.min(U_cg_all[0,99,plta:pltb,plta:pltb]),np.max(U_cg_all[0,99,plta:pltb,plta:pltb])), fontsize=10)
#
    ax2 = plt.subplot(212, sharex=ax1)
    ax2.set_aspect(depth_plot_ratio)
    plt.contourf(XC[icy,plta:pltb]*1e-3,RC.squeeze(),U_theta[:,icy,plta:pltb],vt_range,cmap=cm.Blues_r)
    cb=plt.colorbar(format='%.3f', ticks=vt_ticks)
    cb.set_label(r'$U_\theta \, (m/s)$', labelpad=-40, y=1.1, rotation=0)
    cb.ax.tick_params(labelsize=10)
#    cb.formatter.set_scientific(True)
    CS1 = plt.contour(XC[icy,plta:pltb]*1e-3,RC.squeeze(),U_theta[:,icy,plta:pltb], levels, colors='0.6')
    plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)

#    plt.text(375,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
#    plt.text(375,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title(r'Tangential velocity section, $U_\theta$, timestep %d hr' % (it_hour),fontsize=11)
#
    plt.tight_layout (pad = 1)
    if it_hour==0:
        plt.savefig('./figures/Utheta_0000%d.png' % (it_hour))
    elif it_hour<10:
        plt.savefig('./figures/Utheta_000%d.png' % (it_hour))
    elif it_hour<100:
        plt.savefig('./figures/Utheta_00%d.png' % (it_hour))
    elif it_hour<1000:
        plt.savefig('./figures/Utheta_0%d.png' % (it_hour))
    else:
        plt.savefig('./figures/Utheta_%d.png' % (it_hour))
    plt.close()

##########################################################
    fig = plt.figure(figsize=(5.5,9))
#
    ax1 = plt.subplot(211)
    ax1.set_aspect(1)
    plt.contourf(xc_dom,yc_dom,U_r[idepth,plta:pltb,plta:pltb],vr_range,cmap=cm.seismic)#
    rmax = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
    #axc = plt.gca()
    ax2.add_artist(rmax)
    cb=plt.colorbar(format='%1.3f', ticks=vr_ticks)
    cb.set_label(r'$U_r \, (m/s)$', labelpad=-40, y=1.1, rotation=0)
    cl = plt.contour(xc_dom,yc_dom,U_r[idepth,plta:pltb,plta:pltb], levels, colors='0.6')
    plt.clabel(cl, fmt='%2.2f', colors='k', fontsize=10)
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(r'Radial velocity, $U_r$, timestep %d hr' % (it_hour),fontsize=11)
#
    ax2 = plt.subplot(212, sharex=ax1)
    ax2.set_aspect(depth_plot_ratio)
    plt.contourf(XC[icy,plta:pltb]*1e-3,RC.squeeze(),U_r[:,icy,plta:pltb],vr_range,cmap=cm.seismic)
    cb=plt.colorbar(format='%.3f', ticks=vr_ticks)
    cb.set_label(r'$U_r \, (m/s)$', labelpad=-40, y=1.1, rotation=0)
    cb.ax.tick_params(labelsize=10)
#    cb.formatter.set_scientific(True)
    CS1 = plt.contour(XC[icy,plta:pltb]*1e-3,RC.squeeze(),U_r[:,icy,plta:pltb], levels, colors='0.6')
    plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)

#    plt.text(375,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
#    plt.text(375,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title(r'Radial velocity section, $U_r$, timestep %d hr' % (it_hour),fontsize=11)
#
    plt.tight_layout (pad = 1)
    if it_hour==0:
        plt.savefig('./figures/Ur_0000%d.png' % (it_hour))
    elif it_hour<10:
        plt.savefig('./figures/Ur_000%d.png' % (it_hour))
    elif it_hour<100:
        plt.savefig('./figures/Ur_00%d.png' % (it_hour))
    elif it_hour<1000:
        plt.savefig('./figures/Ur_0%d.png' % (it_hour))
    else:
        plt.savefig('./figures/Ur_%d.png' % (it_hour))
    plt.close()
#

#======= time loop end
#####################################################################
