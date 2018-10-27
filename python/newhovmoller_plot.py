"""
Plot the hovmoller diagrams
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors

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

dumpfreq = 7200
ts = 60  # time step = ts*dt (in second); = 7200 = dumpfreq
endtime = 21600
nit = int(endtime/ts) # no of time step 

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')
dXC = mit.rdmds('DXC')
dYC = mit.rdmds('DYC')

nr = RC.size  # no of grid vertical 
nx = XC[0,:].size # no of grid in x
ny = YC[0,:].size # no of grid in y 

sec_y = int(ny/2)

#Uall=np.ones((nr,ny,nx,nit))
#Vall=np.ones((nr,ny,nx,nit))
Wall=np.ones((nr,ny,nx,nit))

#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(0,nit,1)
dstart = (itrs[0]+1)*ts/720 #5
dend = (itrs[-1]+1)*ts/720  #10
#itrs= [3600]

# time loop
for it in itrs:
#    U = mit.rdmds('U',((it+1)*ts))
    """
    # Sometimes for some reason 1 time step file cannot be read, 
    # this just a "quick solution" to skip the broken timestep
    if (it+1)*30==1080:
        W = mit.rdmds('W',((it+2)*ts))
    else:
        W = mit.rdmds('W',((it+1)*ts))
    """
#    Uall[:,:,:,(it-itrs[0])]= U;
#    U = []
#
#    V = mit.rdmds('V',((it+1)*ts))
#    Vall[:,:,:,(it-itrs[0])]= V;
#    V = []
#
    W = mit.rdmds('W',((it+1)*ts))
    Wall[:,:,:,(it-itrs[0])]= W;
    W = []

plot_perturb = 0
if plot_perturb == 1:
    u_tilde = 0.0*np.ones((nr,ny,nx,nit));
    v_tilde = 0.0*np.ones((nr,ny,nx,nit));
    dyU = 0.0*np.ones((ny,nx,nit));
    dxV = 0.0*np.ones((ny,nx,nit));
    zeta = 0.0*np.ones((ny,nx,nit));

    for it in itrs:
        u_tilde[:,:,:,it] = Uall[:,:,:,it]-Uall[:,:,:,0]
        v_tilde[:,:,:,it] = Vall[:,:,:,it]-Vall[:,:,:,0]
#
        for i in range(1,250):
            for j in range(1,250):
                dyU[j,i,it] = (u_tilde[1,j,i,it]-u_tilde[1,j-1,i,it])/dYC[j,i];
                dxV[j,i,it] = (v_tilde[1,j,i,it]-v_tilde[1,j,i-1,it])/dXC[j,i];
                zeta[j,i,it] = dxV[j,i,it]-dyU[j,i,it];

# Compute Kinetic Energy of perturbation
"""
K=u_tilde**2+v_tilde**2

lnsqrtK=np.log(np.sqrt(K))
"""
#######################Plot the hovmoller diagrams##########################################

## Hovmoller time vs y 
#
time = itrs*dumpfreq/3600
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
#w_range = np.linspace(-0.6,0.6,61, endpoint=True)
#w_range2 = np.linspace(-0.3,0.3,61, endpoint=True)
w_range = 101

hstart1 = 0 
hend1 = 60
"""
plt.figure(figsize=(12,6))
plt.plot(time,lnsqrtK[0,125,125,:])
plt.ylabel(r'$\ln K^{1/2}$')
plt.xlabel("time (hour)")
plt.title(r'$\ln K^{1/2}$ at $x = 150km, y = 150km$, day %d-%d' % (dstart, dend))
plt.savefig('./figures/K.png')

fig = plt.figure(figsize=(15,6))
#
ax1 = fig.add_subplot(1, 2, 1)
plt.contourf(time[hstart1:hend1].squeeze()*2,YC[:,125]*1e-3,Wall[33,:,125,hstart1:hend1]*1e3,w_range,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
plt.colorbar(label='$W \ [mm/s]$', format='%1.3f')
CS1 = plt.contour(time[hstart1:hend1].squeeze()*2,YC[:,125]*1e-3,Wall[33,:,125,hstart1:hend1]*1e3, levels, colors='0.6')
plt.clabel(CS1, fmt='%2.3f', colors='k', fontsize=10)
plt.ylabel("y (km)")
plt.xlabel("time (hour)")
plt.title(r'$W$ at $x = 150km, z\sim 207m$, day %d-%d' % (int(hstart1/12), int(hend1/12)))
#

ax2 = fig.add_subplot(1, 2, 2)
hstart2 = 120
hend2 = 240
plt.contourf(time[hstart2:hend2].squeeze(),YC[:,125]*1e-3,Wall[33,:,125,hstart2:hend2]*1e3,101,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
plt.colorbar(label='$W \ [mm/s]$', format='%1.3f')
CS3 = plt.contour(time[hstart2:hend2].squeeze(),YC[:,125]*1e-3,Wall[33,:,125,hstart2:hend2]*1e3, levels, colors='0.6')
plt.clabel(CS3, fmt='%2.3f', colors='k', fontsize=10)
plt.ylabel("y (km)")
plt.xlabel("time (hour)")
plt.title(r'$W$ at $x = 150km, z\sim 207m$, day %d-%d' % (int(hstart2/24), int(hend2/24)))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_W207_day%d_%d.png' % (dstart, dend))

#### Hovmoller surface  vorticity

fig = plt.figure(figsize=(12,6))
#
plt.contourf(time.squeeze(),YC[:,125]*1e-3,zeta[:,125,:]*1e3,101,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
plt.colorbar(label='$W \ [mm/s]$', format='%1.3f')

CS2 = plt.contour(time.squeeze(),YC[:,125]*1e-3,zeta[:,125,:]*1e3, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)

plt.ylabel("y (km)")
plt.xlabel("time (hour)")
plt.title('Surface vorticity ($\zeta$) at $x=150km$, day %d-%d' % (dstart, dend))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_Vorticity_day%d_%d.png' % (dstart, dend))

#### Hovmoller surface u tilde
fig = plt.figure(figsize=(12,6))
#
plt.contourf(time.squeeze(),YC[:,125]*1e-3,u_tilde[0,:,125,:]*1e3,101,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
plt.colorbar(label='$W \ [mm/s]$', format='%1.3f')

CS2 = plt.contour(time.squeeze(),YC[:,125]*1e-3,u_tilde[0,:,125,:]*1e3, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)

plt.ylabel("y (km)")
plt.xlabel("time (hour)")
plt.title('U perturbation ($\zeta$) at $x=150km$, day %d-%d' % (dstart, dend))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_utilde_day%d_%d.png' % (dstart, dend))
"""
#### Hovmoller time vs depth
fig = plt.figure(figsize=(9,5))
#
plt.contourf(time.squeeze(),RC.squeeze(),Wall[:,sec_y-100,sec_y,:]*1e3,101,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
plt.colorbar(label='$W \ [mm/s]$', format='%1.2f')

CS2 = plt.contour(time.squeeze(),RC.squeeze(),Wall[:,sec_y-100,sec_y,:]*1e3, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)

plt.ylabel("depth (m)")
plt.xlabel("time (hour)")
#plt.title('Vertical velocity ($W$) at $x=%dkm,\  y=%dkm$, day %d-%d' % (XC[sec_y,sec_y]*1e-3, YC[sec_y-20,sec_y]*1e-3,dstart, dend))
plt.title('Vertical velocity ($W$) 24km south of eddy core, day %d-%d' % (dstart, dend))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_Wcoredown_day%d_%d.png' % (dstart, dend))

#
#####################################################################

fig = plt.figure(figsize=(9,5))
#
plt.contourf(time.squeeze(),YC[plta+100:pltb-100,sec_y],Wall[idepth,plta+100:pltb-100,sec_y+20]*1e3,101,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
plt.colorbar(label='$W \ [mm/s]$', format='%1.2f')

CS2 = plt.contour(time.squeeze(),YC[plta+100:pltb-100,sec_y],Wall[idepth,plta+100:pltb-100,sec_y]*1e3, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)

plt.ylabel("y (km)")
plt.xlabel("time (hour)")
#plt.title('Vertical velocity ($W$) at $x=%dkm,\  y=%dkm$, day %d-%d' % (XC[sec_y,sec_y]*1e-3, YC[sec_y-20,sec_y]*1e-3,dstart, dend))
plt.title('Vertical velocity ($W$) 24km south of eddy core, day %d-%d' % (dstart, dend))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_trapped_day%d_%d.png' % (dstart, dend))

