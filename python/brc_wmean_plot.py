"""
Plot the w on the surface and the cross section
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

plt.ion()

timestep = 120 #timestep input in second
dumpfreq = 7200 # file every 3 hours for dumpfreq 10800
file_dit = dumpfreq/timestep
day_s = 10
day_e = 15
startfile = day_s*file_dit #or any integer*file_dit
endfile = day_e*86400/timestep + 50 #1 day 
itrs = np.arange(startfile,endfile, file_dit)
#itrs= [2400,4800,6240] #2400, 3600, 4800,7200]#
nit = itrs.size
XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')

#v_range = 100
v_range_max = 6.5
v_range = np.linspace(-v_range_max, v_range_max, 101, endpoint=True) #.02 or .12 
v_ticks = np.linspace(-v_range_max, v_range_max, 11, endpoint=True)
#levels = np.linspace(-0.02, 0.02,6, endpoint=True)
levels = np.concatenate((np.linspace(-10.0,0,10,endpoint=False),np.linspace(1.0,10.0,10,endpoint=True)),axis=0)
conv_factor = 86400
conv_unit = 'm/day'

idepth = 25 #21 14

plot_depth_ratio = 0.5# 0.7 0.15
Rmax = 25*1e3 #25 15.625*1e3
posext = -520 
posextx = 400 #440 500
nx = XC.shape[1]
ny = YC.shape[0]
nr = RC.size
icx = int(nx/2) 
icy = int(ny/2)

plta = 40
pltb = 290
pltdepth = 65
xc_dom =XC[plta:pltb, plta:pltb]*1e-3
yc_dom =YC[plta:pltb, plta:pltb]*1e-3
nu=0.01 #from kpp=0.01 or Az=4*1e-3
f0=6.5e-5
de = -np.sqrt(2*nu/f0)
XY = np.sqrt(XC[icy,:]**2+YC[:,icx]**2)
Wdiag = np.zeros((nr,nx))
W5d=np.ones((nit,nr,ny,nx))

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

# time loop
for it in itrs:
    W = mit.rdmds('W',it)
    W5d[int((it-itrs[0])/file_dit), :,:,:]= W;

#edit here
#    Wmint = np.trapz(W[: ,: ,:],x=-RC.squeeze(), axis = 0)/(-RC[idepth])
#    for i in np.arange(0,nx):
#        Wdiag[:,i] = W[:,i,i]

W_ta = np.mean(W5d, axis=0)
wmin = np.min(W_ta)*conv_factor
wmax = np.max(W_ta)*conv_factor
#
#==============================================================================================
# plot vertically averaged W and section

fig = plt.figure(figsize=(5.5,9))
ax1 = plt.subplot(211)
ax1.set_aspect(1)
plt.contourf(xc_dom,yc_dom,W_ta[idepth,plta:pltb,plta:pltb]*conv_factor,v_range,cmap=cm.seismic)#
plt.plot(XC[plta:pltb,icx]*1e-3,YC[plta:pltb,icx]*1e-3, ':')#axvline(x=XC[166,166])
cb=plt.colorbar(format='%1.1f', ticks=v_ticks)
cb.set_label(r'$W \, $['+conv_unit+']', labelpad=-40, y=1.1, rotation=0)
cb.ax.tick_params(labelsize=10)
CS2 = plt.contour(xc_dom,yc_dom,W_ta[idepth,plta:pltb,plta:pltb]*conv_factor, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
ax1.add_artist(rmax1)
rmax2 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax1.add_artist(rmax2)
#    plt.text(470,470,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)

plt.xlabel("x (km)", fontsize=10)
plt.ylabel("y (km)", fontsize=10)
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
#    plt.title(r'Vertically averaged $W_{num}$, timestep %d hr' % (int(it*timestep/3600)),fontsize=11)
plt.title(r'$W_{num}(z=%d)$, time averaged day %d to %d' % (RC[idepth], day_s, day_e), fontsize=11)


#
ax2 = plt.subplot(212, sharex=ax1)
ax2.set_aspect(plot_depth_ratio)
plt.contourf(YC[plta:pltb,icx]*1e-3,RC[:pltdepth].squeeze(),W_ta[:pltdepth,plta:pltb,icx]*conv_factor,v_range,cmap=cm.seismic)
cb=plt.colorbar(format='%.1f', ticks=v_ticks)
cb.set_label(r'$W \, $['+conv_unit+']', labelpad=-40, y=1.1, rotation=0)
plt.plot(YC[plta:pltb,icx]*1e-3, np.tile(RC[idepth,0,0], (YC[plta:pltb,icx].size)), ls='dashed', color='0.6')#axvline(x=XC[166,166])
#    plt.plot(YC[plta:pltb,icx]*1e-3, np.tile(de*3, (YC[plta:pltb,icx].size)), ls='dashed', color='0.6')#axvline(x=XC[166,166])
plt.text(560,-80,'$z=%d m$' % (RC[idepth]), fontsize=10, color='0.6')
#    plt.text(580,5,'$3 \delta_E =%d m$' % (de*3), fontsize=10, color='0.6')
cb.ax.tick_params(labelsize=10)
#    cb.formatter.set_scientific(True)
CS1 = plt.contour(YC[plta:pltb,int(icx/2)]*1e-3,RC[:pltdepth].squeeze(),W_ta[:pltdepth,plta:pltb,int(icx/2)]*conv_factor, levels, colors='0.6')
plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)

#    plt.text(375,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
#    plt.text(375,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

plt.xlabel("y (km)")
plt.ylabel("z (m)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
plt.title('$W_{num}$ section, time averaged day %d to %d' % (day_s, day_e), fontsize=11)
plt.text(posextx,posext,'$W_{extrema}=$[%1.1f, %1.1f] ' % (wmin, wmax) + conv_unit, fontsize=12)
plt.tight_layout (pad = 1)
plt.savefig('./figures/W/W_timeav_mday_day%d-%d.png' % (day_s, day_e))
#plt.close()


