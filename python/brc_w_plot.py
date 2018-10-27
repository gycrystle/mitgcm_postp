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
day_s = 0
day_e = 15
startfile = day_s*86400/timestep #or any integer*file_dit
endfile = day_e*86400/timestep + 50 #1 day 
itrs = np.arange(startfile,endfile, file_dit)
#itrs= [10800] #2400, 3600, 4800,7200]#


XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')

#v_range = 100
v_range_max = 1.0
v_range = np.linspace(-v_range_max, v_range_max, 101, endpoint=True) #.02 or .12 
v_ticks = np.linspace(-v_range_max, v_range_max, 11, endpoint=True)
#levels = np.linspace(-0.02, 0.02,6, endpoint=True)
#levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
levels = np.concatenate((np.linspace(-20.0,0,10,endpoint=False),np.linspace(2,20.0,10,endpoint=True)),axis=0)

idepth = 25 #21 14

plot_depth_ratio = 0.5# 0.7 0.15
Rmax = 25*1e3 #25 15.625*1e3
posext = -520 
posextx = 400 #440 500
conv_factor = 86400 # to m/day
conv_unit = 'm/day'
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
#edit here
    Wmint = np.trapz(W[: ,: ,:],x=-RC.squeeze(), axis = 0)/(-RC[idepth])
    wmax=np.max(W)*conv_factor
    wmin=np.min(W)*conv_factor
    for i in np.arange(0,nx):
        Wdiag[:,i] = W[:,i,i]

#
    """
    fig = plt.figure(figsize=(7.5,6))
    plt.contourf(xc_dom,yc_dom,W[idepth,plta:pltb,plta:pltb]*1e3,v_range,cmap=cm.seismic)#
    #plt.plot(XC[:,166]*1e-3,YC[:,166]*1e-3, ':')#axvline(x=XC[166,166])
    plt.colorbar(label="W (mm/s)")
    CS2 = plt.contour(XC[sec_x,plta:pltb]*1e-3,YC[plta:pltb,sec_y]*1e-3,W[idepth,plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
    plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)

    plt.text(480,480,'$W_{extrema}=$[ %1.3f, %1.3f] $mm/s$' % (wmin, wmax))
#    plt.text(375,375,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title('W at $z \sim %d m$, timestep %d' % (RC[idepth],it),fontsize=11)
    if it==0:
        plt.savefig('./figures/W_%d_0000%d.png' % (RC[idepth],it))
    elif it<100:
        plt.savefig('./figures/W_%d_000%d.png' % (RC[idepth],it))
    elif it<1000:
        plt.savefig('./figures/W_%d_00%d.png' % (RC[idepth],it))
    elif it<10000:
        plt.savefig('./figures/W_%d_0%d.png' % (RC[idepth],it))
    else:
        plt.savefig('./figures/W_%d_%d.png' % (RC[idepth],it))
    plt.close()


#
    fig = plt.figure(figsize=(7.5,6))
#
#    ax1 = fig.add_subplot(1, 2, 1)
#    ax1.set_aspect(1)
    plt.contourf(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3,v_range,cmap=cm.seismic)#
    #plt.plot(XC[:,166]*1e-3,YC[:,166]*1e-3, ':')#axvline(x=XC[166,166])
    plt.colorbar(label="W (m/s)")
    CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb], levels, colors='0.6')
    plt.clabel(CS2, fmt='%.2e', colors='k', fontsize=10)

    plt.text(375,390,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
    plt.text(375,375,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d' % (RC[idepth],it),fontsize=11)
    if it==0:
        plt.savefig('./figures/Wi_%d_0000%d.png' % (RC[idepth],it))
    elif it<100:
        plt.savefig('./figures/Wi_%d_000%d.png' % (RC[idepth],it))
    elif it<1000:
        plt.savefig('./figures/Wi_%d_00%d.png' % (RC[idepth],it))
    elif it<10000:
        plt.savefig('./figures/Wi_%d_0%d.png' % (RC[idepth],it))
    else:
        plt.savefig('./figures/Wi_%d_%d.png' % (RC[idepth],it))
#    plt.close()
#
    fig = plt.figure(figsize=(10,6))
#
    #ax1 = fig.add_subplot(1, 2, 1)
    plt.contourf(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3,v_range,cmap=cm.seismic)
    plt.colorbar(label="W (mm/s)")
    CS1 = plt.contour(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3, levels, colors='0.6')
    plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)

    plt.text(375,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
    plt.text(375,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("y (km)")
    plt.ylabel("z (m)")
    plt.title('Vertical velocity $W$ at $x=%dkm$, timestep %d' % (XC[int(si_x/2),int(si_x/2)]*1e-3,it), fontsize=12)
    #
    plt.tight_layout(pad=1)
    if it==0:
        plt.savefig('./figures/W_section_0000%d.png' % (it))
    elif it<100:
        plt.savefig('./figures/W_section_000%d.png' % (it))
    elif it<1000:
        plt.savefig('./figures/W_section_00%d.png' % (it))
    elif it<10000:
        plt.savefig('./figures/W_section_0%d.png' % (it))
    else:
        plt.savefig('./figures/W_section_%d.png' % (it))
    plt.close()
    """
#==============================================================================================
# plot vertically averaged W and section

    fig = plt.figure(figsize=(5.5,9))
    ax1 = plt.subplot(211)
    ax1.set_aspect(1)
    plt.contourf(xc_dom,yc_dom,W[idepth,plta:pltb,plta:pltb]*conv_factor,v_range,cmap=cm.seismic)#
#    plt.contourf(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*conv_factor,v_range,cmap=cm.seismic)#
    plt.plot(XC[plta:pltb,icx]*1e-3,YC[plta:pltb,icx]*1e-3, ':')#axvline(x=XC[166,166])
    cb=plt.colorbar(format='%1.1f', ticks=v_ticks)
    cb.set_label(r'$W \, $['+conv_unit+']', labelpad=-40, y=1.1, rotation=0)
    cb.ax.tick_params(labelsize=10)
    CS2 = plt.contour(xc_dom,yc_dom,W[idepth,plta:pltb,plta:pltb]*conv_factor, levels, colors='0.6')
#    CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*conv_factor, levels, colors='0.6')
    plt.clabel(CS2, fmt='%2.1f', colors='k', fontsize=10)
    rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
    ax1.add_artist(rmax1)
    rmax2 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
    ax1.add_artist(rmax2)
#    plt.text(470,470,'$W_{extrema}=$[%1.3f, %1.3f]' + conv_unit % (wmin, wmax), fontsize=10)
    #plt.text(375,375,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("x (km)", fontsize=10)
    plt.ylabel("y (km)", fontsize=10)
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
#    plt.title(r'Vertically averaged $W_{num}$, timestep %d hr' % (int(it*timestep/3600)),fontsize=11)
    plt.title(r'$W_{num}(z=%d)$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)

#
    ax2 = plt.subplot(212, sharex=ax1)
    ax2.set_aspect(plot_depth_ratio)
    plt.contourf(YC[plta:pltb,icx]*1e-3,RC[:pltdepth].squeeze(),W[:pltdepth,plta:pltb,icx]*conv_factor,v_range,cmap=cm.seismic)
    cb=plt.colorbar(format='%.1f', ticks=v_ticks)
    cb.set_label(r'$W \, $['+conv_unit+']', labelpad=-40, y=1.1, rotation=0)
    plt.plot(YC[plta:pltb,icx]*1e-3, np.tile(RC[idepth,0,0], (YC[plta:pltb,icx].size)), ls='dashed', color='0.6')#axvline(x=XC[166,166])
#    plt.plot(YC[plta:pltb,icx]*1e-3, np.tile(de*3, (YC[plta:pltb,icx].size)), ls='dashed', color='0.6')#axvline(x=XC[166,166])
    plt.text(560,-90,'$z=%d m$' % (RC[idepth]), fontsize=10, color='0.6')
#    plt.text(580,5,'$3 \delta_E =%d m$' % (de*3), fontsize=10, color='0.6')
    cb.ax.tick_params(labelsize=10)
#    cb.formatter.set_scientific(True)
    CS1 = plt.contour(YC[plta:pltb,int(icx/2)]*1e-3,RC[:pltdepth].squeeze(),W[:pltdepth,plta:pltb,icx]*conv_factor, levels, colors='0.6')
    plt.clabel(CS1, fmt='%2.1f', colors='k', fontsize=10)

#    plt.text(375,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
#    plt.text(375,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("y (km)")
    plt.ylabel("z (m)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
    plt.title('$W_{num}$ section, timestep %d hr' % (int(it*timestep/3600)), fontsize=11)
    plt.text(posextx,posext,'$W_{extrema}=$[%1.1f, %1.1f] $m/day$' % (wmin, wmax), fontsize=12)
    plt.tight_layout (pad = 1)
    if it==0:
        plt.savefig('./figures/W/W_0000%d.png' % (it))
    elif it<100:
        plt.savefig('./figures/W/W_000%d.png' % (it))
    elif it<1000:
        plt.savefig('./figures/W/W_00%d.png' % (it))
    elif it<10000:
        plt.savefig('./figures/W/W_0%d.png' % (it))
    else:
        plt.savefig('./figures/W/W_%d.png' % (it))
    plt.close()

# Landscape version

    fig = plt.figure(figsize=(11,4.5))
    ax1 = plt.subplot(121)
    #ax1 = fig.add_subplot(2, 1, 1)
    ax1.set_aspect(1)
    plt.contourf(xc_dom,yc_dom,W[idepth,plta:pltb,plta:pltb]*conv_factor,v_range,cmap=cm.seismic)#
#    plt.contourf(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3,v_range,cmap=cm.seismic)#
    plt.plot(XC[icy, plta:pltb]*1e-3,YC[plta:pltb,icx]*1e-3, ':')#axvline(x=XC[166,166])
    cb=plt.colorbar(format='%1.2f', ticks=v_ticks)
    cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
    CS2 = plt.contour(xc_dom,yc_dom,W[idepth,plta:pltb,plta:pltb]*conv_factor, levels, colors='0.6')
#    CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
    plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
    rmax22 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
    ax1.add_artist(rmax22)
    rmax21 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
    ax1.add_artist(rmax21)

    ax1.set_aspect(1)

#    plt.text(470,470,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
    #plt.text(375,375,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
#    plt.title(r'Vertically averaged $W_{num}$, timestep %d hr' % (int(it*timestep/3600)),fontsize=11)
    plt.title(r'$W_{num}(z=%d)$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)

#
#    ax2 = plt.subplot(122, sharex=ax1)
    ax2 = plt.subplot(122)

    #ax2 = fig.add_subplot(2, 1, 2)
    ax2.set_aspect(plot_depth_ratio+0.2) #0.7 or 0.32
    plt.contourf(XY[plta:pltb]*1e-3,RC[:pltdepth].squeeze(),Wdiag[:pltdepth,plta:pltb]*conv_factor,v_range,cmap=cm.seismic)
    cb=plt.colorbar(format='%.2f', ticks=v_ticks)
    cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
#    plt.plot(XY[plta:pltb]*1e-3, np.tile(de*3, (YC[plta:pltb,icx].size)), ls='dashed', color='0.6')#axvline(x=XC[166,166])
    plt.plot(XY[plta:pltb]*1e-3, np.tile(RC[idepth,0,0], (YC[plta:pltb,icx].size)), ls='dashed', color='0.6')#axvline(x=XC[166,166])
#    plt.text(580,-50,'$3 \delta_E =%d m$' % (de*3), fontsize=10, color='0.6')
    plt.text(800,-50,'$z=%d m$' % (RC[idepth]), fontsize=10, color='0.6')
    CS1 = plt.contour(XY[plta:pltb]*1e-3,RC[:pltdepth].squeeze(),Wdiag[:pltdepth,plta:pltb]*conv_factor, levels, colors='0.6')
    plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)

    plt.xlabel("xy (km)")
    plt.ylabel("z (m)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
    plt.title('$W_{num}$ section, timestep %d hr' % (int(it*timestep/3600)), fontsize=11)
    plt.text(posextx+200,posext+50,'$W_{extrema}=$[%1.3f, %1.3f] $m/day$' % (wmin, wmax), fontsize=12)
    plt.tight_layout (pad = 1)
    ax1.set_aspect(1)

    if it==0:
        plt.savefig('./figures/W/Wl_0000%d.png' % (it))
    elif it<100:
        plt.savefig('./figures/W/Wl_000%d.png' % (it))
    elif it<1000:
        plt.savefig('./figures/W/Wl_00%d.png' % (it))
    elif it<10000:
        plt.savefig('./figures/W/Wl_0%d.png' % (it))
    else:
        plt.savefig('./figures/W/Wl_%d.png' % (it))
    plt.close()

