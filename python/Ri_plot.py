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
dz = np.diff(-RC, axis=0)
dumpfreq = 7200
timestep = 120
idepth = 74
icx=int(XC[0,:].size/2)
icy=int(YC[0,:].size/2)

temp_range = np.linspace(20, 25, 101, endpoint=True)
rho0 = 999.8
f0 = 6.5e-5
alphaK = 2.0e-4
#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
#itrs=np.arange(1170,10980,180)
itrs= [0,60,360,720, 3600,7200]

outicks = [1e7, 1e11, 1e15, 1e27, 1e31]

#fig, ax1 = plt.subplots(figsize=(7,6))
plt.figure(figsize=(7,6))
ax1 = plt.subplot(121)
plt.xlabel(r'$Ri$')
plt.ylabel("Depth [m]")
plt.title('Inside eddy, %d km from core' % (abs(XC[145,145]-XC[165,165])*1e-3))
plt.grid(True)
ax2 = plt.subplot(122, sharey=ax1)
#logs = ax2.coords[0]
# time loop
for it in itrs:
    Temp = mit.rdmds('T',it)
    rho = rho0*(1-alphaK*Temp)
    N2i = 9.8*np.diff(rho[:,145,145], axis=0)/dz.squeeze()/rho0
    N2o = 9.8*np.diff(rho[:,45,45], axis=0)/dz.squeeze()/rho0
    U = mit.rdmds('U', it)
    dzU2i = (np.diff(U[:,145,145], axis=0)/dz.squeeze())**2
    dzU2o = (np.diff(U[:,45,45], axis=0)/dz.squeeze())**2
    V = mit.rdmds('V', it)
    dzV2i = (np.diff(V[:,145,145], axis=0)/dz.squeeze())**2
    dzV2o = (np.diff(V[:,45,45], axis=0)/dz.squeeze())**2
    Rii = N2i/(dzU2i+dzV2i)
    Rio = N2o/(dzU2o+dzV2o)
#    lnRi = np.log(Ri)
#
#    plt.figure(figsize=(4.5,6))
    ax1.semilogx(Rii[:65], RC[1:66].squeeze())#, label='timestep %d hr' % (it*timestep/3600))
    ax2.semilogx(Rio[:65], RC[1:66].squeeze(), label='timestep %d hr' % (it*timestep/3600))
    #plt.plot(rhop[:,165,165], -zc, label='core Initial Cond')
    #plt.plot(rho_anomd[49:,0,0], -zc[49:], label='@Rmax 3D-cons')
    #plt.plot(rhop[:,165,165], -zc, label='@Rmax Initial Cond')
#
plt.xlabel(r'$Ri$')
plt.grid(True)
plt.ylim(-100,0)
#plt.ylabel("Depth [m]")
plt.title('Outside eddy, %d km from core' % (abs(XC[45,45]-XC[165,165])*1e-3))
plt.tight_layout(pad=1)
ax2.legend(loc='best')
#logs = ax2.coords[0]
#logs.set_ticks(outicks)
#plt.savefig('./figures/Richardson%d.png' %(it))
plt.savefig('./figures/Richardson10day.png')


#    visc = mit.rdmds('KPPviscAz',it)
    
