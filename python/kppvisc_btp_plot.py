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

#ts = 60  # time step = ts*dt (in second); 
timestep = 180 #timestep input in second
dumpfreq = 10800 # file every 3 hours for dumpfreq 10800
file_dit = dumpfreq/timestep

#ts = file_dit
file_dit = 720
day_s = 0
day_e = 10

startfile = day_s*file_dit+180 #or any integer*file_dit
endfile = day_e*86400/timestep + 50 #1 day 
itrs = np.arange(startfile,endfile, file_dit)

ipdepth = 99
idepth = 33 # depth for subsurface drifts 36 for approx -99 in 75 grid 

#itrs = np.arange(0, 21610, file_dit)
itrs = [60,480,2400,4800,7200]#, 14400, 18000, 21600] #, 2160, 2880, 3600, 4320, 5040, 5760, 6480, 7200]

nit = len(itrs) #itrs.size

rho0 = 999.8
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
depth_plot_ratio= 0.8#0.7 or 0.5
plta = 30
pltb = 150
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
         
ke_range = np.linspace(0,0.27,100, endpoint=True)
ke_ticks = np.linspace(0,0.27,11, endpoint=True)
#ke_range = 99
#ke_ticks = 11

plt.figure(figsize=(7,6))
ax1 = plt.subplot(121)
plt.xlabel(r'$A_v$')
plt.ylabel("Depth [m]")
plt.title('Inside eddy, at core')
plt.grid(True)
ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
plt.xlabel(r'$A_v$')
plt.grid(True)
plt.ylim(-300,0)
#plt.ylabel("Depth [m]")
plt.title('Outside eddy, %d km from core' % (abs(XC[50,50]-XC[90,90])*1e-3))

for it in itrs:
    viscKPP = mit.rdmds('KPPviscAz',it)
#    V = mit.rdmds('V',it)
#
#

#    ax1.semilogx(viscKPP[:65,165,165], RC[:65].squeeze())#, label='timestep %d hr' % (it*timestep/3600))
#    ax2.semilogx(viscKPP[:65,70,70], RC[:65].squeeze(), label='timestep %d hr' % (it*timestep/3600))

#############Plot figures ##############
############# visc ####################################
    fig = plt.figure(figsize=(11,4))
#
    ax3 = fig.add_subplot(1, 2, 1)
    ax3.set_aspect(1)
    plt.contourf(xc_d,yc_d,viscKPP[idepth,plta:pltb,plta:pltb],ke_range,cmap=cm.ocean_r)
    plt.colorbar(label=r'$\nu_z  m^2s^{-2}$', format='%.2f', ticks=ke_ticks)
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(r'$\nu_z m^2s^{-1}$,z=%d m, time step %d hr' % (RC[idepth], int(it*timestep/3600)))
#
    ax4 = fig.add_subplot(1, 2, 2)
    ax4.set_aspect(depth_plot_ratio)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf (YC[plta:pltb, sec_y]*1e-3, RC[:ipdepth].squeeze(), viscKPP[:ipdepth,plta:pltb, sec_y], ke_range,cmap=cm.ocean_r)
    plt.colorbar(label=r'$\nu_z \ [m^2s^{-1}]$', format='%.2f', ticks=ke_ticks)
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title(r'$\nu_z m^2s^{-2}$, time step %d hr' % (int(it*timestep/3600)))
#
    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/visckpp_0000"+ str(int(it)) + ".png")
    elif it < 100:
        plt.savefig("./figures/visckpp_000"+ str(int(it)) + ".png")
    elif it < 1000:
        plt.savefig("./figures/visckpp_00"+ str(int(it)) + ".png")
    elif it < 10000:
        plt.savefig("./figures/visckpp_0"+ str(int(it)) + ".png")
    else:
        plt.savefig("./figures/visckpp_"+ str(int(it)) + ".png")
    plt.close()
    """
    plt.figure(figsize=(4.5,6))
    plt.semilogx(viscKPP[:,90,90], RC.squeeze())
    #plt.plot(rhop[:,165,165], -zc, label='core Initial Cond')
    #plt.plot(rho_anomd[49:,0,0], -zc[49:], label='@Rmax 3D-cons')
    #plt.plot(rhop[:,165,165], -zc, label='@Rmax Initial Cond')
    plt.xlabel(r'$\nu_z$ [m^2/s]')
    plt.ylabel("Depth [m]")
    plt.title("KPP viscosity vertical")
    plt.tight_layout(pad=1)
    plt.legend(loc='best')
    plt.savefig('KPP_Az_%d.png' %(it))
    """
#    ax1.semilogx(viscKPP[:65,165,165], RC[:65].squeeze())#, label='timestep %d hr' % (it*timestep/3600))
#    ax2.semilogx(viscKPP[:65,90,90], RC[:65].squeeze(), label='timestep %d hr' % (it*timestep/3600))
    ax1.plot(viscKPP[:,90,90], RC.squeeze())#, label='timestep %d hr' % (it*timestep/3600))
    ax2.plot(viscKPP[:,50,50], RC.squeeze(), label='timestep %d hr' % (it*timestep/3600))


#plt.xlabel(r'$A_v$')
#plt.grid(True)
#plt.ylim(-300,0)
#plt.ylabel("Depth [m]")
#plt.title('Outside eddy, %d km from core' % (abs(XC[90,90]-XC[165,165])*1e-3))
plt.tight_layout(pad=1)
ax2.legend(loc='lower right')
#logs = ax2.coords[0]
#logs.set_ticks(outicks)
#plt.savefig('./figures/Richardson%d.png' %(it))
plt.savefig('./figures/KPPvisc_profile.png')

    
"""
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
