
"""
Python script to compare dissipation of vorticity in various running cases
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors

plt.ion()

def vorticity(U,V,dXC,dYC):
    nx = U[0,0,:].size
    ny = U[0,:,0].size
    nr = U[:,0,0].size
    dyU  = np.zeros((nr, ny,nx));
    dxV  = np.zeros((nr, ny,nx));
    zeta = np.zeros((nr, ny,nx));
    for i in range(1,nx):
        for j in range(1,ny):
            dyU[:,j,i] = (U[:,j,i]-U[:,j-1,i])/dYC[j,i];
            dxV[:,j,i] = (V[:,j,i]-V[:,j,i-1])/dXC[j,i];
            zeta[:,j,i] = dxV[:,j,i]-dyU[:,j,i];
    return zeta

# select plot domain
plta = 40
pltb = 290
idepth = 21
posext = 400
depth_plot_ratio= 0.7# 7

#vort_range = np.linspace(-0.55, 0.55, 101, endpoint=True)
#vorticks = np.linspace(-0.55, 0.55, 11, endpoint=True)
#zero_contour = [0]

XC  = mit.rdmds('XC')
dXC = mit.rdmds('DXC')
YC  = mit.rdmds('YC')
dYC = mit.rdmds('DYC')
RC  = mit.rdmds('RC')


nr = RC.size
nx = XC[0,:].size
ny = YC[:,0].size

timestep = 120 #timestep input in second
dumpfreq = 7200 # in second

day_s = 0
day_e = 5
ts = dumpfreq/timestep  # time step = ts*dt (in second); = 7200 = dumpfreq

startfile = day_s*ts #or any integer*file_dit
endfile = day_e*86400/timestep + 50 #1 day 

itrs = np.arange(startfile,endfile, ts)
nit=itrs.size
#nit=len(itrs)

xc_d = XC[plta:pltb,plta:pltb]*1e-3
yc_d = YC[plta:pltb,plta:pltb]*1e-3
icx = int(nx/2)
icy = int(ny/2)

svortc0 = np.zeros((nit))
svortc1 = np.zeros((nit))
svortc2 = np.zeros((nit))
svortc3 = np.zeros((nit))
svortc4 = np.zeros((nit))
svortc5 = np.zeros((nit))

rho0 = 999.8
f0 = 6.5e-5
f00 = 1e-4

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

time = itrs*timestep/3600

#Plot vorticity profile
itrsd = [0,720,2160,3600,5760, 7200]
fig = plt.figure(figsize=(8,6))
ax1 = plt.subplot(121)
plt.xlabel(r'$\zeta/f_0$')
plt.ylabel("Depth [m]")
plt.title('Vorticity at initial eddy core, case Av = 1e-2')
plt.grid(True)
ax2 = plt.subplot(122,  sharey=ax1)
plt.xlabel(r'$W m/day$')
#plt.ylabel("Depth [m]")
plt.title('W at %d km from initial eddy core' %((XC[145,145]-XC[165,165])*1e-3))
plt.grid(True)


for it in itrsd:
    U = mit.rdmds('U',it)
    V = mit.rdmds('V',it)
    zetaz = vorticity(U, V, dXC, dYC)
    W = mit.rdmds('W',it)

    """
    plt.figure (figsize = (4.5, 6))
    plt.plot(np.tile(1e-5,(RC.size)), RC.squeeze(), label=run0_label)
    plt.plot(np.tile(4e-3,(RC.size)), RC.squeeze(), label=run1_label)
    plt.plot(np.tile(1e-2,(RC.size)), RC.squeeze(), label=run2_label)
    plt.plot(viscKPP3[:,90,90], RC.squeeze(),ls='--', label=run3_label)
    plt.plot(np.tile(4e-3,(RC.size)), RC.squeeze(), label=run4_label)
    plt.plot(viscKPP5[:,90,90], RC200.squeeze(),ls='--', label=run5_label)
#    plt.plot(np.tile(1e-5,(RC.size)), RC.squeeze())
#    plt.plot(viscKPP0[:,90,90], RC.squeeze())
    plt.xlabel(r'$\nu_z [m^2/s]$')
    plt.ylabel("Depth [m]")
    plt.title("viscosity vertical")
    plt.tight_layout(pad=1)
    plt.legend(loc='best')
    plt.savefig('./figures/Compare_Az_%d.png' %(it))
    """
#
    ax1.plot(zetaz[:60,165,165]/f0, RC[:60].squeeze(), label='timestep %d hr' % (it*timestep/3600))
    ax2.plot(W[:60,145,145]*86400, RC[:60].squeeze(), label='timestep %d hr' % (it*timestep/3600))

ax1.legend(loc='lower left')
plt.tight_layout(pad=1)
plt.savefig('./figures/W_Vorticity60.png')

"""
# surface vort core difference
fig = plt.figure(figsize=(7,4))
plt.plot(time,-(svortc1-svortc0)/f0, label='anticyclone')
plt.plot(time,(svortc2+svortc0)/f0, label='cyclone')
plt.xlabel("time (hour)")
plt.ylabel(r'$\frac{|\zeta_c|-|\zeta_c^{ref}|}{|f|}$', fontsize=14)
plt.legend()
plt.xlim(0,itrs[-1]*timestep/3600)
plt.title(r'Wind modified core vorticity, $\frac{|\zeta_c|-|\zeta_c^{ref}|}{|f|}$', fontsize=14)
plt.tight_layout(pad=1)
plt.savefig("./figures/modcore.png")

"""
