"""
Python script to read compute and plot Vorticity (zeta) output from MITgcm
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
#import cv2 as cv

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
def center_of_mass(X):
    # calculate center of mass of a closed polygon
    x = X[:,0]
    y = X[:,1]
    g = (x[:-1]*y[1:] - x[1:]*y[:-1])
    A = 0.5*g.sum()
    cx = ((x[:-1] + x[1:])*g).sum()
    cy = ((y[:-1] + y[1:])*g).sum()
    return 1./(6*A)*np.array([cx,cy])

nit = 60 # no of time step 
ts = 60  # time step = ts*dt (in second); 
timestep = 120 #timestep input in second
dumpfreq = 7200 # file every 3 hours for dumpfreq 10800
file_dit = dumpfreq/timestep
day_s = 0
day_e = 15
startfile = day_s*file_dit #or any integer*file_dit
endfile = day_e*86400/timestep + 50 #1 day 
#itrs = np.arange(startfile,endfile, file_dit)


vort_range = np.linspace(-0.7, 0.7, 101, endpoint=True)
vorticks = np.linspace(-0.7, 0.7, 11, endpoint=True)
zero_contour = [0]
#vort_range = 101

#itrs = np.arange(0, 7210, 2400)
itrs = [0,2400,4800, 7200] #, 2160, 2880, 3600, 4320, 5040, 5760, 6480, 7200]

nit = len(itrs) #itrs.size
rho0 = 999.8
f0 = 1e-4
XC = mit.rdmds('XC')
dXC = mit.rdmds('DXC')
YC = mit.rdmds('YC')
dYC = mit.rdmds('DYC')
RC = mit.rdmds('RC')
XG = mit.rdmds('XG')
YG = mit.rdmds('YG')

nr = RC.size
nx = XC[0,:].size
ny = YC[:,0].size

U_c = np.zeros((nr,ny,nx))
V_c = np.zeros((nr,ny,nx))
svortc = np.zeros((nit))

# values for plots
depth_plot_ratio= 0.15#7
plta = 40
pltb = 290
posext = 400
xc_d = XC[plta:pltb,plta:pltb]*1e-3
yc_d = YC[plta:pltb,plta:pltb]*1e-3
sec_y = int(nx/2)

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

dyU = np.zeros((nr,ny,nx));
dxV = np.zeros((nr,ny,nx));
         
zeta = np.zeros((nr,ny,nx));
ke_range = np.linspace(0,8,100)
ke_range = 99
fig1 = plt.figure(figsize=(4.5,3))
figt = fig1.add_subplot(1, 1, 1)
plt.xlabel("x (km)")
plt.ylabel(r'$\frac{\zeta_0}{|f|}$')
plt.title(r'Surface vorticity section')

for it in itrs:
    U = mit.rdmds('U',it)
    V = mit.rdmds('V',it)
#
    for i in range(1,nx):
        for j in range(1,ny):
            #for k in range(0,nr):
            dyU[:,j,i] = (U[:,j,i]-U[:,j-1,i])/dYC[j,i];
            dxV[:,j,i] = (V[:,j,i]-V[:,j,i-1])/dXC[j,i];
            zeta[:,j,i] = dxV[:,j,i]-dyU[:,j,i];

# Must interpolate to cell center before computing the tangential and radial velocity
    for z in np.arange(0,nr):
        for y in np.arange(0,ny):
            U_c[z,y,:] = np.interp(XC[y,:],XG[y,:],U[z,y,:])
        for x in np.arange(0,nx):
            V_c[z,:,x] = np.interp(YC[:,x],YG[:,x],V[z,:,x])

    KE_m =(U_c**2 + V_c**2)/2
    svortc [int((it-itrs[0])/2400)]=np.min(zeta[0,:,:]) #file_dit
#    KE_m =np.sum((U**2 + V**2), axis=0)
#
    zetamin = np.min(zeta)/f0
#############Plot figures ##############
    figt.plot(XC[sec_y,plta:pltb]*1e-3,zeta[0,sec_y,plta:pltb]/f0, label='t=%d hr' %(int(it*timestep/3600)))
    #plt.figure()
    fig = plt.figure(figsize=(12,4))
#
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    plt.contourf(xc_d,yc_d,zeta[0,plta:pltb,plta:pltb]/f0, vort_range,cmap=cm.seismic)
#    plt.contourf(xc_d,yc_d,zeta[0,plta:pltb,plta:pltb]/f0, vort_range,vmin=zetamin,vmax=-zetamin,cmap=cm.seismic)
    cb=plt.colorbar(format='%.2f', ticks=vorticks)
    cb.set_label(r'$\frac{\zeta}{|f|}$', labelpad=-40, y=1.1, rotation=0)
    CS2 = plt.contour(XC[0,plta+50:pltb-50]*1e-3,YC[plta+50:pltb-50,0]*1e-3,zeta[0,plta+50:pltb-50,plta+50:pltb-50]/f0, zero_contour, colors='0.6')
    c =  center_of_mass(CS2.allsegs[-1][0])
    ax1.plot(c[0],c[1], marker="+", markersize=10, color='0.6')
#    plt.show()
#    M = cv.moments(CS2)
#    cx = int(M['m10']/M['m00'])
#    cy = int(M['m01']/M['m00'])
#    plt.plot(cx,cy, marker="+", markersize=10, color='0.6')
#    plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
#
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(r'Surface vorticity, $\frac{\zeta}{f}$')
#
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_aspect(depth_plot_ratio)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf (XC[sec_y,plta:pltb]*1e-3, RC.squeeze(), zeta[:, sec_y,plta:pltb]/f0, vort_range,cmap=cm.seismic)
    cb=plt.colorbar(format='%.2f', ticks=vorticks)
    cb.set_label(r'$\frac{\zeta}{|f|}$', labelpad=-40, y=1.1, rotation=0)
    plt.text(posext,posext,r'$\frac{\zeta_{min}}{|f|}=$ %1.3f' % (np.min(zeta/f0)))
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("Vorticity")

    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/vort_0000"+ str(it) + ".png")
    elif it < 100:
        plt.savefig("./figures/vort_000"+ str(it) + ".png")
    elif it < 1000:
        plt.savefig("./figures/vort_00"+ str(it) + ".png")
    elif it < 10000:
        plt.savefig("./figures/vort_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/vort_"+ str(it) + ".png")
#    plt.close()

    fig = plt.figure(figsize=(4.7,4))
#
#    ax1 = fig.add_subplot(1, 2, 1)
#    ax1.set_aspect(1)
    plt.contourf(xc_d,yc_d,zeta[0,plta:pltb,plta:pltb]/f0, vort_range,cmap=cm.seismic)
#    plt.contourf(xc_d,yc_d,zeta[0,plta:pltb,plta:pltb]/f0, vort_range,vmin=zetamin,vmax=-zetamin,cmap=cm.seismic)
#    plt.colorbar(label=r'$\frac{\zeta}{|f|}$')
    cb=plt.colorbar(format='%.2f', ticks=vorticks)
    cb.set_label(r'$\frac{\zeta}{|f|}$', labelpad=-40, y=1.1, rotation=0)
    CS2 = plt.contour(XC[0,plta+50:pltb-50]*1e-3,YC[plta+50:pltb-50,0]*1e-3,zeta[0,plta+50:pltb-50,plta+50:pltb-50]/f0, zero_contour, colors='0.6')
    c =  center_of_mass(CS2.allsegs[-1][0])
    plt.plot(c[0],c[1], marker="+", markersize=10, color='0.6')
    plt.text(posext,posext,r'$\frac{\zeta_{min}}{|f|}=$ %1.3f' % (np.min(zeta/f0)))
    plt.show()
#    M = cv.moments(CS2)
#    cx = int(M['m10']/M['m00'])
#    cy = int(M['m01']/M['m00'])
#    plt.plot(cx,cy, marker="+", markersize=10, color='0.6')
#    plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
#
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(r'Surface vorticity, $\frac{\zeta}{|f|}$, time step %d hr' % (int(it*timestep/3600)))
#
    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/svort_0000"+ str(it) + ".png")
    elif it < 100:
        plt.savefig("./figures/svort_000"+ str(it) + ".png")
    elif it < 1000:
        plt.savefig("./figures/svort_00"+ str(it) + ".png")
    elif it < 10000:
        plt.savefig("./figures/svort_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/svort_"+ str(it) + ".png")
#    plt.close()

############# KE ####################################
    fig = plt.figure(figsize=(12,4))
#
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    plt.contourf(xc_d,yc_d,KE_m[0,plta:pltb,plta:pltb],ke_range,cmap=cm.ocean_r)
    plt.colorbar(label=r'$KE  m^2s^{-2}$')
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(r'KE $m^2s^{-2}$')
#
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_aspect(depth_plot_ratio)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf (XC[sec_y,plta:pltb]*1e-3, RC.squeeze(), KE_m[:, sec_y,plta:pltb], ke_range,cmap=cm.ocean_r)
    plt.colorbar(label=r'$KE \ [m^2s^{-2}]$')
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("KE")

    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/KEc_0000"+ str(it) + ".png")
    elif it < 100:
        plt.savefig("./figures/KEc_000"+ str(it) + ".png")
    elif it < 1000:
        plt.savefig("./figures/KEc_00"+ str(it) + ".png")
    elif it < 10000:
        plt.savefig("./figures/KEc_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/KEc_"+ str(it) + ".png")
    plt.close()
# surface vort core
plt.tight_layout(pad=1)
fig = plt.figure(figsize=(5,3))
plt.plot(itrs*timestep/3600,svortc/f0)
plt.xlabel("time (hour)")
plt.ylabel(r'$\frac{\zeta}{|f|}$')
plt.title(r'Minimum Surface vorticity, $\frac{\zeta_{min}}{|f|}$')
plt.tight_layout(pad=1)
plt.savefig("./figures/Surfvortmin.png")


##========== from diag ====================================
#itrsd = np.arange(7200, 7210, 720)
itrsd = [0]
for it in itrsd:
#    vort = mit.rdmds('momvort',it)
    KE = mit.rdmds('momKE',it)
    """
    fig = plt.figure(figsize=(12,4))    

    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf(xc_d,yc_d,vort[0,plta:pltb,plta:pltb]*1e5, vort_range,MidpointNormalize(vmax=-vmin, midpoint=0.0),cmap=cm.seismic)
    plt.colorbar(label=r'$\zeta \ 10^5 s^{-1}$')
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(r'Surface vorticity, $\zeta \ 10^5 s^{-1}$')
    #plt.savefig("mld"+ str(i) + ".png")

    ax7 = fig.add_subplot(1, 2, 2)
    ax7.set_aspect(depth_plot_ratio)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf(XC[sec_y,plta:pltb]*1e-3,RC.squeeze(),vort[:,sec_y,plta:pltb]/f0, vort_range,cmap=cm.seismic)
    plt.colorbar(label=r'$\zeta \ 10^5 s^{-1}$')
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("Vorticity")

    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/vortmom_000"+ str(it) + ".png")
    elif it < 1000:
        plt.savefig("./figures/vortmom_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/vortmom_"+ str(it) + ".png")
#    plt.close()
    """

# KE
    """
    fig = plt.figure(figsize=(12,4))    

    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf(xc_d,yc_d,KE[0,plta:pltb,plta:pltb], ke_range,cmap=cm.ocean_r)
    plt.colorbar(label='KE')
    plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title("KE")
    #plt.savefig("mld"+ str(i) + ".png")

    ax7 = fig.add_subplot(1, 2, 2)
    ax7.set_aspect(0.7)
    #plt.subplot(121, aspect='equal', adjustable='box-forced')
    plt.contourf(XC[sec_y,plta:pltb]*1e-3,RC.squeeze(),KE[:,sec_y,plta:pltb], ke_range,cmap=cm.ocean_r)
    plt.colorbar(label=r'$\zeta$')
    plt.xlabel("x (km)")
    plt.ylabel("z (m)")
    plt.title("KE")

    plt.tight_layout(pad=1)
    if it ==0:
        plt.savefig("./figures/KE_000"+ str(it) + ".png")
    elif it < 1000:
        plt.savefig("./figures/KE_0"+ str(it) + ".png")
    else:
        plt.savefig("./figures/KE_"+ str(it) + ".png")
#    plt.close()
    """
"""
#plt.contour(xi, yi, zi, v, linewidths=0.5, colors='k')
#plt.contourf(xi, yi, zi, v, cmap=plt.cm.jet)
#x = plt.colorbar(ticks=v)
"""
