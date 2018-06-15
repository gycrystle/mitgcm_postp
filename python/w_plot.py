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
day_e = 5
startfile = day_s*file_dit #or any integer*file_dit
endfile = day_e*86400/timestep + 100 #1 day 
#itrs = np.arange(startfile,endfile, file_dit)
itrs= [0,60,120]#

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')

#v_range = 100
v_range = np.linspace(-0.08, 0.08, 101, endpoint=True)
v_ticks = np.linspace(-0.08, 0.08, 11, endpoint=True)
#levels = np.linspace(-0.02, 0.02,6, endpoint=True)
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)

idepth = 24
si_x = XC.shape[1]
si_y = YC.shape[0]
sec_x = int(si_x/2)
icx =int(si_x/2) 
sec_y = int(si_y/2)
plta = 40
pltb = 290
xc_dom =XC[plta:pltb, plta:pltb]*1e-3
yc_dom =YC[plta:pltb, plta:pltb]*1e-3

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
    wmax=np.max(W)*1e3
    wmin=np.min(W)*1e3
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
#    plt.close()


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
#    plt.close()
    """
#==============================================================================================
# plot vertically averaged W and section

    fig = plt.figure(figsize=(6,9))
    ax1 = plt.subplot(211)
    #ax1 = fig.add_subplot(2, 1, 1)
    ax1.set_aspect(1)
    plt.contourf(xc_dom,yc_dom,W[idepth,plta:pltb,plta:pltb]*1e3,v_range,cmap=cm.seismic)#
#    plt.contourf(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3,v_range,cmap=cm.seismic)#
    plt.plot(XC[plta:pltb,icx]*1e-3,YC[plta:pltb,icx]*1e-3, ':')#axvline(x=XC[166,166])
    cb=plt.colorbar(format='%1.3f', ticks=v_ticks)
    cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
    cb.ax.tick_params(labelsize=10)
    CS2 = plt.contour(xc_dom,yc_dom,W[idepth,plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#    CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
    plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
#    plt.text(470,470,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
    #plt.text(375,375,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("x (km)", fontsize=10)
    plt.ylabel("y (km)", fontsize=10)
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
#    plt.title(r'Vertically averaged $W_{num}$, timestep %d hr' % (int(it*timestep/3600)),fontsize=11)
    plt.title(r'$W_{num}(z=%d)$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)

#
    ax2 = plt.subplot(212, sharex=ax1)
    #ax2 = fig.add_subplot(2, 1, 2)
    ax2.set_aspect(0.7)
    plt.contourf(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)],v_range,cmap=cm.seismic)
    cb=plt.colorbar(format='%.3f', ticks=v_ticks)
    cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
    cb.ax.tick_params(labelsize=10)
#    cb.formatter.set_scientific(True)
    CS1 = plt.contour(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3, levels, colors='0.6')
    plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)

#    plt.text(375,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
#    plt.text(375,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("y (km)")
    plt.ylabel("z (m)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
    plt.title('$W_{num}$ section, timestep %d hr' % (int(it*timestep/3600)), fontsize=11)
    plt.text(470,-320,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
    plt.tight_layout (pad = 1)
    if it==0:
        plt.savefig('./figures/W_0000%d.png' % (it))
    elif it<100:
        plt.savefig('./figures/W_000%d.png' % (it))
    elif it<1000:
        plt.savefig('./figures/W_00%d.png' % (it))
    elif it<10000:
        plt.savefig('./figures/W_0%d.png' % (it))
    else:
        plt.savefig('./figures/W_%d.png' % (it))
#    plt.close()

# Landscape version

    fig = plt.figure(figsize=(11,4.8))
    ax1 = plt.subplot(121)
    #ax1 = fig.add_subplot(2, 1, 1)
    ax1.set_aspect(1)
    plt.contourf(xc_dom,yc_dom,W[idepth,plta:pltb,plta:pltb]*1e3,v_range,cmap=cm.seismic)#
#    plt.contourf(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3,v_range,cmap=cm.seismic)#
    plt.plot(XC[plta:pltb,icx]*1e-3,YC[plta:pltb,icx]*1e-3, ':')#axvline(x=XC[166,166])
    cb=plt.colorbar(format='%1.3f', ticks=v_ticks)
    cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
    CS2 = plt.contour(xc_dom,yc_dom,W[idepth,plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#    CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
    plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
    ax1.set_aspect(1)

#    plt.text(470,470,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
    #plt.text(375,375,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
#    plt.title(r'Vertically averaged $W_{num}$, timestep %d hr' % (int(it*timestep/3600)),fontsize=11)
    plt.title(r'$W_{num}(z=%d)$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)

#
    ax2 = plt.subplot(122, sharex=ax1)
    #ax2 = fig.add_subplot(2, 1, 2)
    ax2.set_aspect(0.7) #0.7 or 0.32
    plt.contourf(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3,v_range,cmap=cm.seismic)
    cb=plt.colorbar(format='%.3f', ticks=v_ticks)
    cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
    CS1 = plt.contour(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3, levels, colors='0.6')
    plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)

#    plt.text(375,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
#    plt.text(375,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("y (km)")
    plt.ylabel("z (m)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
    plt.title('$W_{num}$ section, timestep %d hr' % (int(it*timestep/3600)), fontsize=11)
    plt.text(470,-320,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
    plt.tight_layout (pad = 1)
    ax1.set_aspect(1)

    if it==0:
        plt.savefig('./figures/Wl_0000%d.png' % (it))
    elif it<100:
        plt.savefig('./figures/Wl_000%d.png' % (it))
    elif it<1000:
        plt.savefig('./figures/Wl_00%d.png' % (it))
    elif it<10000:
        plt.savefig('./figures/Wl_0%d.png' % (it))
    else:
        plt.savefig('./figures/Wl_%d.png' % (it))
#    plt.close()

