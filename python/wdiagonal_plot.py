"""
Plot the w on the surface and the cross section
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm

#plt.ion()

timestep = 120 #timestep input in second
dumpfreq = 7200 # file every 3 hours for dumpfreq 10800
file_dit = dumpfreq/timestep
day_s = 0
day_e = 15
startfile = day_s*file_dit #or any integer*file_dit
endfile = day_e*86400/timestep + 100 #1 day 
itrs = np.arange(startfile,endfile, file_dit)
#itrs= [10800]#

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')
XY = np.sqrt(XC[165,:]**2+YC[:,165]**2)

#v_range = 100
v_range = np.linspace(-0.05, 0.05, 101, endpoint=True)
#levels = np.linspace(-0.02, 0.02,6, endpoint=True)
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)

idepth = 24
si_x = XC.shape[1]
si_y = YC.shape[0]
sec_x = int(si_x/2)
sec_y = int(si_y/2)
plta = 10
pltb = 320
xc_dom =XC[plta:pltb, plta:pltb]*1e-3
yc_dom =YC[plta:pltb, plta:pltb]*1e-3
XY = np.sqrt(XC[sec_x,:]**2+YC[:,sec_x]**2)
Wdiag = np.zeros((nr,nx))

# time loop
for it in itrs:
    W = mit.rdmds('W',it)
#edit here
    Wmint = np.trapz(W[: ,: ,:],x=-RC.squeeze(), axis = 0)/(-RC[idepth])
    wmax=np.max(W)*1e3
    wmin=np.min(W)*1e3
    for i in np.arange(0,nx):
        Wdiag[:,i] = W[:,i,i]
#
#==============================================================================================
# plot vertically averaged W and section

    fig = plt.figure(figsize=(6,9))
    ax1 = plt.subplot(211)
    #ax1 = fig.add_subplot(2, 1, 1)
    #ax1.set_aspect(1)
    plt.contourf(xc_dom,yc_dom,W[idepth,plta:pltb,plta:pltb]*1e3,v_range,cmap=cm.seismic)#
#    plt.contourf(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3,v_range,cmap=cm.seismic)#
    plt.plot(XC[plta:pltb,165]*1e-3,YC[plta:pltb,165]*1e-3, ':')#axvline(x=XC[166,166])
    cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
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
    ax2 = plt.subplot(212, sharex=ax1)
    #ax2 = fig.add_subplot(2, 1, 2)
    ax2.set_aspect(0.3)
    plt.contourf(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3,v_range,cmap=cm.seismic)
    cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
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
    plt.text(350,-1550,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
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
    plt.close()

# Landscape version + diagonal cross section
    fig = plt.figure(figsize=(11.5,5))
    ax1 = plt.subplot(121)
    ax1.set_aspect(1)
    plt.contourf(xc_dom,yc_dom,W[idepth,plta:pltb,plta:pltb]*1e3,v_range,cmap=cm.seismic)#
#    plt.contourf(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3,v_range,cmap=cm.seismic)#
#    plt.plot(XC[plta:pltb,165]*1e-3,YC[plta:pltb,165]*1e-3, ':')#axvline(x=XC[166,166])
    plt.plot(XC[165,plta:pltb]*1e-3,YC[plta:pltb,165]*1e-3, ':')#axvline(x=XC[166,166])
    cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
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
    """
    ax2 = plt.subplot(132)
    #ax2 = fig.add_subplot(2, 1, 2)
    ax2.set_aspect(0.3) #0.7 or 0.32
    plt.contourf(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3,v_range,cmap=cm.seismic)
#    plt.contourf(XY[plta:pltb]*1e-3,RC.squeeze(),Wdiag[:,plta:pltb]*1e3,v_range,cmap=cm.seismic)
#    cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
#    cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
    CS1 = plt.contour(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3, levels, colors='0.6')
    plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)

#    plt.text(375,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
#    plt.text(375,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("y (km)")
    plt.ylabel("z (m)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
    plt.title('$W_{num}$ section, timestep %d hr' % (int(it*timestep/3600)), fontsize=11)
    plt.text(350,-1550,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
    """
###########################
    ax3 = plt.subplot(122)
    ax3.set_aspect(0.5) #0.7 or 0.32
#    plt.contourf(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3,v_range,cmap=cm.seismic)
    plt.contourf(XY[plta:pltb]*1e-3,RC.squeeze(),Wdiag[:,plta:pltb]*1e3,v_range,cmap=cm.seismic)
#    cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
#    cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
    CS1 = plt.contour(XY[plta:pltb]*1e-3,RC.squeeze(),Wdiag[:,plta:pltb]*1e3, levels, colors='0.6')
    plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)

#    plt.text(375,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
#    plt.text(375,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

    plt.xlabel("xy (km)")
    plt.ylabel("z (m)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
    plt.title('$W_{num}$ diagonal section, timestep %d hr' % (int(it*timestep/3600)), fontsize=11)
    plt.text(320,-1380,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)

    plt.tight_layout (pad = 1)
    ax1.set_aspect(1)

    if it==0:
        plt.savefig('./figures/Wd_0000%d.png' % (it))
    elif it<100:
        plt.savefig('./figures/Wd_000%d.png' % (it))
    elif it<1000:
        plt.savefig('./figures/Wd_00%d.png' % (it))
    elif it<10000:
        plt.savefig('./figures/Wd_0%d.png' % (it))
    else:
        plt.savefig('./figures/Wd_%d.png' % (it))
    plt.close()

