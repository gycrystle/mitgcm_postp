"""
Python script to read compute and plot Vorticity (zeta) output from MITgcm
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
#import cv2 as cv

plt.ion()

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)
    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))
def center_of_mass(X):
    # calculate center of mass of a closed polygon
    x = X[:,0]
    y = X[:,1]
    g = (x[:-1]*y[1:] - x[1:]*y[:-1])
    A = 0.5*g.sum()
    cx = ((x[:-1] + x[1:])*g).sum()
    cy = ((y[:-1] + y[1:])*g).sum()
    return 1./(6*A)*np.array([cx,cy])

#nit = 60 # no of time step 
#ts = 60  # time step = ts*dt (in second); 
timestep = 120 #timestep input in second
dumpfreq = 7200 # file every 3 hours for dumpfreq 10800
file_dit = dumpfreq/timestep

#ts = file_dit
#file_dit = 3600
day_s = 0
day_e = 15

startfile = day_s*file_dit #or any integer*file_dit
endfile = day_e*86400/timestep + 50 #1 day 
itrs = np.arange(startfile,endfile, file_dit)

ipdepth = 60
idepth = 25 # depth for subsurface drifts 36 for approx -99 in 75 grid 

#itrs = np.arange(0, 21610, file_dit)
#itrs = [0, 2400,4800, 7200] #, 2160, 2880, 3600, 4320, 5040, 5760, 6480, 7200]

nit = itrs.size
cx = np.zeros((nit))
cy = np.zeros((nit))
c2x = np.zeros((nit))
c2y = np.zeros((nit))

rho0 = 999.8
f0 = 1e-4
XC = mit.rdmds('XC')
dXC = mit.rdmds('DXC')
YC = mit.rdmds('YC')
dYC = mit.rdmds('DYC')
RC = mit.rdmds('RC')
XG = mit.rdmds('XG')
YG = mit.rdmds('YG')

nr = RC.size
nx = XC[0,:].size
ny = YC[:,0].size

U_c = np.zeros((nr,ny,nx))
V_c = np.zeros((nr,ny,nx))
KE_tot = np.zeros((nit))

# values for plots
depth_plot_ratio= 0.5#0.7 or 0.5
plta = 40
pltb = 290
xc_d = XC[plta:pltb,plta:pltb]*1e-3
yc_d = YC[plta:pltb,plta:pltb]*1e-3
sec_y = int(nx/2)
x_c=XC[sec_y,sec_y]*1e-3
y_c=YC[sec_y,sec_y]*1e-3
#pltf = plta#110
pltc = 120#100
pltd = 210#210
plte = pltb#220
xc_dc = XC[pltc:pltd,pltc:pltd]*1e-3
yc_dc = YC[pltc:pltd,pltc:pltd]*1e-3


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


         
ke_range = np.linspace(0,0.1,100, endpoint=True)
ke_ticks = np.linspace(0,0.1,11, endpoint=True)
#ke_range = 99

for it in itrs:
    U = mit.rdmds('U',it)
    V = mit.rdmds('V',it)
#
# Must interpolate to cell center before computing the tangential and radial velocity
    for z in np.arange(0,nr):
        for y in np.arange(0,ny):
            U_c[z,y,:] = np.interp(XC[y,:],XG[y,:],U[z,y,:])
        for x in np.arange(0,nx):
            V_c[z,:,x] = np.interp(YC[:,x],YG[:,x],V[z,:,x])

    KE = (U_c**2 + V_c**2)#/2
    KE_tot[int((it-itrs[0])/file_dit)] = np.sum(KE)
#



#############Plot figures ##############
############# KE ####################################
    fig = plt.figure(figsize=(11,4))
#
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    plt.contourf(xc_d,yc_d,KE[0,plta:pltb,plta:pltb],ke_range,cmap=cm.ocean_r)
    plt.colorbar(label=r'$KE  m^2s^{-2}$', format='%.2f', ticks=ke_ticks)
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(r'surface KE $m^2s^{-2}$, time step %d hr' % (int(it*timestep/3600)))
#
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_aspect(depth_plot_ratio)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf (XC[sec_y,plta:pltb]*1e-3, RC[:ipdepth].squeeze(), KE[:ipdepth, sec_y,plta:pltb], ke_range,cmap=cm.ocean_r)
    plt.colorbar(label=r'$KE \ [m^2s^{-2}]$', format='%.2f', ticks=ke_ticks)
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title(r'KE $m^2s^{-2}$, time step %d hr' % (int(it*timestep/3600)))
#
    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/KE_0000"+ str(int(it)) + ".png")
    elif it < 100:
        plt.savefig("./figures/KE_000"+ str(int(it)) + ".png")
    elif it < 1000:
        plt.savefig("./figures/KE_00"+ str(int(it)) + ".png")
    elif it < 10000:
        plt.savefig("./figures/KE_0"+ str(int(it)) + ".png")
    else:
        plt.savefig("./figures/KE_"+ str(int(it)) + ".png")
    plt.close()

    fig = plt.figure(figsize=(5,4))
#
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.set_aspect(1)
    plt.contourf(xc_d,yc_d,KE[idepth,plta:pltb,plta:pltb],ke_range,cmap=cm.ocean_r)
    plt.colorbar(label=r'$KE  m^2s^{-2}$', format='%.2f', ticks=ke_ticks)
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(r'KE (at z=%d m) $m^2s^{-2}$, time step %d hr' % (RC[idepth], int(it*timestep/3600)))
#
    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/KE50_0000"+ str(int(it)) + ".png")
    elif it < 100:
        plt.savefig("./figures/KE50_000"+ str(int(it)) + ".png")
    elif it < 1000:
        plt.savefig("./figures/KE50_00"+ str(int(it)) + ".png")
    elif it < 10000:
        plt.savefig("./figures/KE50_0"+ str(int(it)) + ".png")
    else:
        plt.savefig("./figures/KE50_"+ str(int(it)) + ".png")
    plt.close()
    

time = itrs*timestep/3600


# KE_evolution
#fig = plt.figure(figsize=(10,5))
fig, ax1 = plt.subplots(figsize=(10,5))
ax1.plot(time,KE_tot)
#ax1.plot(time,cy-cy[0], label='meridional drift')
plt.xlabel("time (hour)")
plt.ylabel('KE')
plt.xlim(0,itrs[-1]*timestep/3600)
#plt.legend(loc='lower left')
plt.title(r'Total kinetic energy evolution')
plt.tight_layout(pad=1)
#ax1.set_aspect(1)
plt.savefig('./figures/time_evol_KE%d.png' % (day_e))

"""
fig, ax1 = plt.subplots(figsize=(5,6))
#plt.plot(time,cx, label='zonal drift')
ax1.plot(cx-cx[0],cy-cy[0], label='surface')
ax1.plot(c2x-c2x[0],c2y-c2y[0], label='%d m' % (RC[idepth]))
plt.xlabel("x (km)")
plt.ylabel('y (km)')
plt.legend(loc='best')
plt.xlim(-4,4)
plt.ylim(-11,0)
ax1.set_aspect(1)
plt.text(cx[-1]-cx[0],cy[-1]-cy[0],'%.2f , %.2f' % (cx[-1]-cx[0],cy[-1]-cy[0]))
plt.text(c2x[-1]-c2x[0],c2y[-1]-c2y[0],'%.2f , %.2f' % (c2x[-1]-c2x[0],c2y[-1]-c2y[0]))
plt.title(r'Eddy drift from center due to wind forcing in %d days' % (day_e-day_s))
plt.tight_layout(pad=1)
plt.savefig('./figures/coretraject%d.png' % (day_e))



# surface vort core
plt.tight_layout(pad=1)
fig = plt.figure(figsize=(5,3))
plt.plot(itrs*timestep/3600,svortc/f0)
plt.xlabel("time (hour)")
plt.ylabel(r'$\frac{\zeta}{|f|}$')
plt.xlim(0,itrs[-1]*timestep/3600)
plt.title(r'Minimum Surface vorticity, $\frac{\zeta_{min}}{|f|}$')
plt.tight_layout(pad=1)
plt.savefig('./figures/Surfvortmin%d.png' % (day_e))

"""
