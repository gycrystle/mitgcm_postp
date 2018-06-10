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
day_e = 15

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
KE_theta = np.zeros((nit,nr,ny,nx))
KE_r = np.zeros((nit,nr,ny,nx))

xc_dom =XC[plta:pltb, plta:pltb]*1e-3
yc_dom =YC[plta:pltb, plta:pltb]*1e-3

dstart = day_s #5
dend = day_e # (itrs[-1]+1)*ts/720  #10
time = itrs*dumpfreq/3600

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
    KE_theta_it = 0.5*U_theta**2
    KE_r_it = 0.5*U_r**2
    
    KE_theta[(it-itrs[0]),:,:,:]= KE_theta_it #
    KE_r[(it-itrs[0]),:,:,:]= KE_r_it #
#======= time loop end

KE_theta_tot = np.trapz(KE_theta,-RC.squeeze(), axis=1)
KE_theta_tot = np.trapz(KE_theta_tot,YC[:,icx], axis=1)
KE_theta_tot = np.trapz(KE_theta_tot,XC[icy,:], axis=1)
#
KE_r_tot = np.trapz(KE_r,-RC.squeeze(), axis=1)
KE_r_tot = np.trapz(KE_r_tot,YC[:,icx], axis=1)
KE_r_tot = np.trapz(KE_r_tot,XC[icy,:], axis=1)
#
#    it_hour = int(it*dumpfreq/3600)
##### plot figure of KE evolution
# Kinetic Energy 
plt.figure(figsize=(12,6))
plt.plot(time,KE_theta_tot, label=r'$KE_\theta$')
#plt.plot(time,KEag[:,int(nx/2+20),int(nx/2+20)], label='ageostrophic')
plt.plot(time,KE_r_tot, label=r'$KE_r$')
plt.plot(time,(KE_theta_tot + KE_r_tot), label=r'$KE_total$')
plt.ylabel(r'$\propto KE$')
plt.xlabel("time (hour)")
plt.legend(loc='best', fontsize=10)
plt.title(r'$\propto KE$ total, day %d-%d' % (dstart, dend))
#    fig = plt.figure(figsize=(5.5,9))
#
#    ax1 = plt.subplot(211)
#    ax1.set_aspect(1)
plt.xlabel("time (hour)")
plt.ylabel(r'$KE (m^2s^{-2})$')
#    plt.title(r'Kinetic Energy' )#% (RC[idepth], it_hour),fontsize=11)

plt.tight_layout (pad = 1)
plt.savefig('./figures/KE_time_evol.png')
plt.close()


#####################################################################
