
"""
Python script to read compute and plot Vorticity (zeta) output from MITgcm
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors

plt.ion()

# select plot domain
plta = 30
pltb = 150
idepth = 99
posext = 500
depth_plot_ratio= 0.7#7

vort_range = np.linspace(-0.7, 0.7, 101, endpoint=True)
vorticks = np.linspace(-0.7, 0.7, 11, endpoint=True)
zero_contour = [0]

XC = mit.rdmds('./run21_00btp/XC')
dXC = mit.rdmds('./run21_00btp/DXC')
YC = mit.rdmds('./run21_00btp/YC')
dYC = mit.rdmds('./run21_00btp/DYC')
RC = mit.rdmds('./run21_00btp/RC')

nr = RC.size
nx = XC[0,:].size
ny = YC[:,0].size

timestep = 180 #timestep input in second
dumpfreq = 10800 # file every 3 hours for dumpfreq 10800

day_s = 0
day_e = 15
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

dyU0 = np.zeros((nr,ny,nx));
dxV0 = np.zeros((nr,ny,nx));
zeta0 = np.zeros((nr,ny,nx));
svortc0 = np.zeros((nit))

dyU1 = np.zeros((nr,ny,nx));
dxV1 = np.zeros((nr,ny,nx));
zeta1 = np.zeros((nr,ny,nx));
svortc1 = np.zeros((nit))

dyU2 = np.zeros((nr,ny,nx));
dxV2 = np.zeros((nr,ny,nx));
zeta2 = np.zeros((nr,ny,nx));
svortc2 = np.zeros((nit))

rho0 = 999.8
f0 = 1e-4

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

#Plot surface vorticity
"""
fig1 = plt.figure(figsize=(4.5,3))
figt = fig1.add_subplot(1, 1, 1)
plt.xlabel("x (km)")
plt.ylabel(r'$\frac{\zeta_0}{|f|}$')
plt.title(r'Surface vorticity section')
"""
for it in itrs:
    U0 = mit.rdmds('./run21_00btp/U',it)
    V0 = mit.rdmds('./run21_00btp/V',it)
    U1 = mit.rdmds('./run22_10btp_kpp/U',it)
    V1 = mit.rdmds('./run22_10btp_kpp/V',it)
    U2 = mit.rdmds('./run25_10btp_cycl/U',it)
    V2 = mit.rdmds('./run25_10btp_cycl/V',it)
#
    for i in range(1,nx):
        for j in range(1,ny):
            #for k in range(0,nr):
            dyU0[:,j,i] = (U0[:,j,i]-U0[:,j-1,i])/dYC[j,i];
            dxV0[:,j,i] = (V0[:,j,i]-V0[:,j,i-1])/dXC[j,i];
            zeta0[:,j,i] = dxV0[:,j,i]-dyU0[:,j,i];
#
            dyU1[:,j,i] = (U1[:,j,i]-U1[:,j-1,i])/dYC[j,i];
            dxV1[:,j,i] = (V1[:,j,i]-V1[:,j,i-1])/dXC[j,i];
            zeta1[:,j,i] = dxV1[:,j,i]-dyU1[:,j,i];
#
            dyU2[:,j,i] = (U2[:,j,i]-U2[:,j-1,i])/dYC[j,i];
            dxV2[:,j,i] = (V2[:,j,i]-V2[:,j,i-1])/dXC[j,i];
            zeta2[:,j,i] = dxV2[:,j,i]-dyU2[:,j,i];

    svortc0 [int(it/ts)]=np.min(zeta0[0,:,:]) #file_dit
    core0=np.where(zeta0==np.min(zeta0[0,:,:]))
    svortc1 [int(it/ts)]=np.min(zeta1[0,:,:]) #file_dit
    svortc2 [int(it/ts)]=np.max(zeta2[0,:,:]) #file_dit
    core2=np.where(zeta2[0,:,:]==np.min(zeta2[0,:,:]))

#    zetamin = np.min(zeta)/f0    
#    figt.plot(XC[sec_y,plta:pltb]*1e-3,zeta[0,sec_y,plta:pltb]/f0, label='t=%d hr' %(int(it*timestep/3600)))
time = itrs*timestep/3600

# surface vort core
fig = plt.figure(figsize=(5,3))
plt.plot(time,-svortc0/f0, label='reference')
plt.plot(time,-svortc1/f0, label='anticyclone')
plt.plot(time,svortc2/f0, label='cyclone')
plt.xlabel("time (hour)")
plt.ylabel(r'$\frac{|\zeta_c|}{|f|}$')
plt.legend()
plt.xlim(0,itrs[-1]*timestep/3600)
plt.title(r'Surface vorticity @core , $\frac{\zeta_{c}}{|f|}$')
plt.tight_layout(pad=1)
plt.savefig("./figures/Surfvortmin.png")

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


