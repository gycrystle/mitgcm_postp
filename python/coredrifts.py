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
 
ts = 60  # time step = ts*dt (in second); 
timestep = 180 #timestep input in second
dumpfreq = 10800 # file every 3 hours for dumpfreq 10800
file_dit = dumpfreq/timestep
#file_dit = 480
day_s = 0
day_e = 15
startfile = day_s*file_dit #or any integer*file_dit
endfile = day_e*86400/timestep + 50 #1 day 
itrs = np.arange(startfile,endfile, file_dit)
ipdepth = 99
idepth = 33 # depth for subsurface drifts 33 or  36 for approx -99 in 75 grid 
vort_range = np.linspace(-0.7, 0.7, 101, endpoint=True)
vorticks = np.linspace(-0.7, 0.7, 11, endpoint=True)
#zero_contour = [0]
#vort_range = 101

#itrs = np.arange(0, 10810, file_dit)
#itrs = [0, 2400,4800, 7200] #, 2160, 2880, 3600, 4320, 5040, 5760, 6480, 7200]

nit = itrs.size
cx = np.zeros((nit))
cy = np.zeros((nit))
c2x = np.zeros((nit))
c2y = np.zeros((nit))

rho0 = 999.8
f0 = 6.5e-5
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
depth_plot_ratio= 0.7#0.7 or 0.5
plta = 30
pltb = 150
xc_d = XC[plta:pltb,plta:pltb]*1e-3
yc_d = YC[plta:pltb,plta:pltb]*1e-3
sec_y = int(nx/2)
x_c=XC[sec_y,sec_y]*1e-3
y_c=YC[sec_y,sec_y]*1e-3
pltf = plta # 0#110
pltc = plta # 30#100
pltd = pltb #150#210
plte = pltb #150#220
xc_dc = XC[pltf:pltd,pltc:pltd]*1e-3
yc_dc = YC[pltf:pltd,pltc:pltd]*1e-3


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

    svortc[int((it-itrs[0])/file_dit)]=np.min(zeta[0,:,:]) #file_dit
#    KE_m =np.sum((U**2 + V**2), axis=0)
#
    zetamin = np.min(zeta)/f0
    zetamax = np.max(zeta)/f0
    if (abs(zetamin)>zetamax):
        zero_contour = [zetamin*3/5]
    else:
        zero_contour = [zetamax*3/5]
#############Plot figures ##############
#    figt.plot(XC[sec_y,plta:pltb]*1e-3,zeta[0,sec_y,plta:pltb]/f0, label='t=%d hr' %(int(it*timestep/3600)))
    #plt.figure()
    fig = plt.figure(figsize=(12,4))
#
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_aspect(1)
    plt.contourf(xc_d,yc_d,zeta[0,plta:pltb,plta:pltb]/f0, vort_range,cmap=cm.seismic)
#    plt.contourf(xc_d,yc_d,zeta[0,plta:pltb,plta:pltb]/f0, vort_range,vmin=zetamin,vmax=-zetamin,cmap=cm.seismic)
    cb=plt.colorbar(format='%.2f', ticks=vorticks)
    cb.set_label(r'$\frac{\zeta}{|f|}$', labelpad=-40, y=1.1, rotation=0)
    CS2 = plt.contour(xc_dc,yc_dc,zeta[0,pltf:pltd,pltc:pltd]/f0, zero_contour, colors='0.6')
    CS4 = plt.contour(xc_dc,yc_dc,zeta[idepth,pltf:pltd,pltc:pltd]/f0, zero_contour, colors='1')
    c =  center_of_mass(CS2.allsegs[0][0])
    cx[int(it/file_dit)]=c[0]
    cy[int(it/file_dit)]=c[1]
    c2 =  center_of_mass(CS4.allsegs[0][0])
    c2x[int(it/file_dit)]=c2[0]
    c2y[int(it/file_dit)]=c2[1]
    ax1.plot(c[0],c[1], marker="+", markersize=8, color='0.6')
#
    plt.close()
#
time = itrs*timestep/3600
figt.legend()
fig1.tight_layout(pad=1)
#fig1.savefig("./figures/vort_section.png")

# core drift
#fig = plt.figure(figsize=(10,5))
fig, ax1 = plt.subplots(figsize=(10,5))
ax1.plot(time,cx-cx[0], label='zonal drift')
ax1.plot(time,cy-cy[0], label='meridional drift')
plt.xlabel("time (hour)")
plt.ylabel('drift in km')
plt.xlim(0,itrs[-1]*timestep/3600)
plt.legend(loc='lower left')
plt.title(r'Eddy drift due to wind forcing')
plt.tight_layout(pad=1)
#ax1.set_aspect(1)
plt.tight_layout(pad=1)
plt.savefig("./figures/coredrift.png")

fig, ax1 = plt.subplots(figsize=(5,6))
#plt.plot(time,cx, label='zonal drift')
ax1.plot(cx-cx[0],cy-cy[0], label='surface')
ax1.plot(c2x-c2x[0],c2y-c2y[0], label='%d m' % (RC[idepth]))
#ax1.arrow(-0.0,0.0,0.0,-24.6, ls='dashed')# label=r'$\dfrac{M_y}{H}=24.6 \ km$')
#ax1.annotate("", xy=(0.0, -24.6), xytext=(0, -25),arrowprops=dict(arrowstyle="<-"))
plt.xlabel("x (km)")
plt.ylabel('y (km)')
plt.legend(loc='best')
plt.xlim(-4,4)
plt.ylim(-27,0)
ax1.set_aspect(0.4)
plt.text(cx[-1]-cx[0],cy[-1]-cy[0],'%.2f , %.2f' % (cx[-1]-cx[0],cy[-1]-cy[0]))
plt.text(c2x[-1]-c2x[0],c2y[-1]-c2y[0],'%.2f , %.2f' % (c2x[-1]-c2x[0],c2y[-1]-c2y[0]))
#plt.text(-3.8,-4.1,r'$30 days * M_y/H$')
#plt.text(-3.8,-5.2,r'$=24.6 \ km$')
plt.title(r'Eddy drift from center due to wind forcing in %d days' % (day_e-day_s))
plt.tight_layout(pad=1)
plt.savefig("./figures/coretraject3_5.png")



# surface vort core
plt.tight_layout(pad=1)
fig = plt.figure(figsize=(5,3))
plt.plot(itrs*timestep/3600,svortc/f0)
plt.xlabel("time (hour)")
plt.ylabel(r'$\frac{\zeta}{|f|}$')
plt.xlim(0,itrs[-1]*timestep/3600)
plt.title(r'Minimum Surface vorticity, $\frac{\zeta_{min}}{|f|}$')
plt.tight_layout(pad=1)
plt.savefig("./figures/Surfvortmin.png")


