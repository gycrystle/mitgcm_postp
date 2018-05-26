"""
Plot the w on the surface and the cross section
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm

plt.ion()

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')

v_range = 100
#v_range = np.linspace(-0.1, 0.1, 91, endpoint=True)
#levels = np.linspace(-0.02, 0.02,6, endpoint=True)
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
#itrs=np.arange(0,10810,60)
itrs= [120]#3600, 10800]#,180,360,720,780,900,1440]#[0,60,120,180,240,300,360,420,480,540,600,660,720]
idepth = 50
si_x = XC.shape[1]
si_y = YC.shape[0]
sec_x = int(si_x/2)
sec_y = int(si_y/2)
plta = 0
pltb = -1
xc_dom =XC[plta:pltb, plta:pltb]*1e-3
yc_dom =YC[plta:pltb, plta:pltb]*1e-3

# time loop
for it in itrs:
    W = mit.rdmds('W',it)
#edit here
    Wmint = np.trapz(W[: ,: ,:],x=-RC.squeeze(), axis = 0)/(-RC[idepth])
    wmax=np.max(W)*1e3
    wmin=np.min(W)*1e3
#
    fig = plt.figure(figsize=(7.5,6))
    plt.contourf(xc_dom,yc_dom,W[idepth,plta:pltb,plta:pltb]*1e3,v_range,cmap=cm.seismic)#
    #plt.plot(XC[:,166]*1e-3,YC[:,166]*1e-3, ':')#axvline(x=XC[166,166])
    plt.colorbar(label="W (mm/s)")
    CS2 = plt.contour(XC[sec_x,plta:pltb]*1e-3,YC[plta:pltb,sec_y]*1e-3,W[idepth,plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
    plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)

    plt.text(375,390,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
    plt.text(375,375,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

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
    plt.colorbar(label="W (mm/s)")
    CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
    plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)

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

