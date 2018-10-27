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

res = 2 # 1 for btp 2 for brc

if res==1:
    dumpfreq = 10800 # 7200 or 10800
    endtime = 7200 # 10800 or 7200

    plta = 30 # 40 or30
    pltb= 150 # 150 or 290

    idepth = 17 # 21 or 14
    ideptht = 66 # depth for trapped 49 fo -199.5m or 66
    ipdepth = 99 # 65 or 99
    iout = 20 # 100 or 48
    iin = 2 # 20 or 10
else:
    dumpfreq = 7200 # 7200 or 10800
    endtime = 10800 # 10800 or 7200

    plta = 40 # 40 or30
    pltb= 290 # 150 or 290

    idepth = 25 # 21 or 14
    ideptht = 49 # depth for trapped 49 fo -199.5m or 66
    ipdepth = 65 # 65 or 99
    iout = 41 # 100 or 48
    iin = 4 # 20 or 10

wa = 12.0 # 0.15 for mm/s
we = 2.0 #
conv_factor = 86400

ts = 60  # time step = ts*dt (in second); = 7200 = dumpfreq
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

Wall=np.ones((nr,ny,nx,nit))

#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(0,nit,1)
dstart = (itrs[0]+1)*dumpfreq/86400 #5
dend = (itrs[-1]+1)*dumpfreq/86400  #10
#itrs= [3600]

time = itrs*dumpfreq/3600
levels = np.concatenate((np.linspace(-20.0,0,10,endpoint=False),np.linspace(2.0,20.0,10,endpoint=True)),axis=0)
we_range = np.linspace(-we,we,101, endpoint=True)
we_ticks = np.linspace(-we,we,11, endpoint=True)
wa_range = np.linspace(-wa,wa,101, endpoint=True)
wa_ticks = np.linspace(-wa,wa,11, endpoint=True)

#w_range2 = np.linspace(-0.3,0.3,61, endpoint=True)
w_range = 101

# time loop
for it in itrs:
#
    W = mit.rdmds('W',((it+1)*ts))
    Wall[:,:,:,(it-itrs[0])]= W;
    W = []

#######################Plot the hovmoller diagrams##########################################

#### Hovmoller time vs depth

fig = plt.figure(figsize=(9,5))
#
plt.contourf(time.squeeze(),RC[:ipdepth].squeeze(),Wall[:ipdepth,sec_y+iin,sec_y,:]*conv_factor,wa_range,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
cb=plt.colorbar(format='%.1f', ticks=wa_ticks)
cb.set_label(r'$W \, [m/day]$', labelpad=-40, y=1.1, rotation=0)

CS2 = plt.contour(time.squeeze(),RC[:ipdepth].squeeze(),Wall[:ipdepth,sec_y+iin,sec_y,:]*conv_factor, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.1f', colors='k', fontsize=10)

plt.ylabel("depth (m)")
plt.xlabel("time (hour)")
plt.title('Vertical velocity ($W$) %dkm north of center of domain, day %d-%d' % ((YC[sec_y+iin, sec_y+iin]-YC[sec_y,sec_y])*1e-3, dstart, dend))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_Wcoreup5_mday_day%d_%d.png' % (dstart, dend))
## ============================================
fig = plt.figure(figsize=(9,5))
#
plt.contourf(time.squeeze(),RC[:ipdepth].squeeze(),Wall[:ipdepth,sec_y-7*iin,sec_y,:]*conv_factor,wa_range,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
cb=plt.colorbar(format='%.1f', ticks=wa_ticks)
cb.set_label(r'$W \, [m/day]$', labelpad=-40, y=1.1, rotation=0)

CS2 = plt.contour(time.squeeze(),RC[:ipdepth].squeeze(),Wall[:ipdepth,sec_y-3*iin,sec_y,:]*conv_factor, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.1f', colors='k', fontsize=10)

plt.ylabel("depth (m)")
plt.xlabel("time (hour)")
plt.title('Vertical velocity ($W$) %dkm south of center of domain, day %d-%d' % ((YC[sec_y-7*iin, sec_y-7*iin]-YC[sec_y,sec_y])*1e-3, dstart, dend))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_Wcoredown5_mday_day%d_%d.png' % (dstart, dend))
## ===============================================
fig = plt.figure(figsize=(9,5))
#
plt.contourf(time.squeeze(),RC[:ipdepth].squeeze(),Wall[:ipdepth,sec_y+iout,sec_y,:]*conv_factor,we_range,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
cb=plt.colorbar(format='%.1f', ticks=we_ticks)
cb.set_label(r'$W \, [m/day]$', labelpad=-40, y=1.1, rotation=0)

CS2 = plt.contour(time.squeeze(),RC[:ipdepth].squeeze(),Wall[:ipdepth,sec_y+iout,sec_y,:]*conv_factor, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.1f', colors='k', fontsize=10)

plt.ylabel("depth (m)")
plt.xlabel("time (hour)")
#plt.title('Vertical velocity ($W$) at $x=%dkm,\  y=%dkm$, day %d-%d' % (XC[sec_y,sec_y]*1e-3, YC[sec_y-20,sec_y]*1e-3,dstart, dend))
plt.title('Vertical velocity ($W$) %dkm north of center of domain, day %d-%d' % ((YC[sec_y+iout,sec_y+iout]-YC[sec_y,sec_y])*1e-3,dstart, dend))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_Wout5_mday_day%d_%d.png' % (dstart, dend))
#plt.close()

#
#####################################################################
## Hovmoller time vs y 

fig = plt.figure(figsize=(9,5))
#
plt.contourf(time.squeeze(),YC[plta:pltb,sec_y]*1e-3,Wall[idepth,plta:pltb,sec_y]*conv_factor,wa_range,cmap=cm.seismic)
cb=plt.colorbar(format='%.1f', ticks=wa_ticks)
cb.set_label(r'$W \, [m/day]$', labelpad=-40, y=1.1, rotation=0)

CS2 = plt.contour(time.squeeze(),YC[plta:pltb,sec_y]*1e-3,Wall[idepth,plta:pltb,sec_y]*conv_factor, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.1f', colors='k', fontsize=10)

plt.ylabel("y (km)")
plt.xlabel("time (hour)")
#plt.title('Vertical velocity ($W$) at $x=%dkm,\  y=%dkm$, day %d-%d' % (XC[sec_y,sec_y]*1e-3, YC[sec_y-20,sec_y]*1e-3,dstart, dend))
plt.title('Vertical velocity ($W$) at z=%dm, day %d-%d' % (RC[idepth], dstart, dend))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_eklayer_mday_day%d_%d.png' % (dstart, dend))

fig = plt.figure(figsize=(9,5))
#
plt.contourf(time.squeeze(),YC[plta:pltb,sec_y]*1e-3,Wall[ideptht,plta:pltb,sec_y]*conv_factor,wa_range,cmap=cm.seismic)
cb=plt.colorbar(format='%.1f', ticks=wa_ticks)
cb.set_label(r'$W \, [m/day]$', labelpad=-40, y=1.1, rotation=0)

CS2 = plt.contour(time.squeeze(),YC[plta:pltb,sec_y]*1e-3,Wall[ideptht,plta:pltb,sec_y]*conv_factor, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.1f', colors='k', fontsize=10)

plt.ylabel("y (km)")
plt.xlabel("time (hour)")
#plt.title('Vertical velocity ($W$) at $x=%dkm,\  y=%dkm$, day %d-%d' % (XC[sec_y,sec_y]*1e-3, YC[sec_y-20,sec_y]*1e-3,dstart, dend))
plt.title('Vertical velocity ($W$) at z=%dm, day %d-%d' % (RC[ideptht], dstart, dend))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_z200_mday_day%d_%d.png' % (dstart, dend))


