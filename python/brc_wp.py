"""
Compute 
W from pressure gradient

"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.ion()

# ====================================================================
# ========= Functions and Formatting =================================
def center_of_mass(X):
    # calculate center of mass of a closed polygon
    x = X[:,0]
    y = X[:,1]
    g = (x[:-1]*y[1:] - x[1:]*y[:-1])
    A = 0.5*g.sum()
    cx = ((x[:-1] + x[1:])*g).sum()
    cy = ((y[:-1] + y[1:])*g).sum()
    return 1./(6*A)*np.array([cx,cy])

def vorticity(U,V,dXC,dYC):
    nx = U[0,:].size
    ny = U[:,0].size
    dyU = np.zeros((ny,nx));
    dxV = np.zeros((ny,nx));
    zeta = np.zeros((ny,nx));
    for i in range(1,nx):
        for j in range(1,ny):
            dyU[j,i] = (U[j,i]-U[j-1,i])/dYC[j,i];
            dxV[j,i] = (V[j,i]-V[j,i-1])/dXC[j,i];
            zeta[j,i] = dxV[j,i]-dyU[j,i];
    return zeta

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

# =======================================================================    

# select plot domain
plta = 40
pltb = 290
idepth = 65 #depth to substract geostrophic flow 
posext = 500 # location of text 
posexty = 400
timestep = 120 #timestep input in second
dumpfreq = 7200 # in second

day_s = 10
day_e = 15

ts = dumpfreq/timestep  # time step = ts*dt (in second); = 7200 = dumpfreq
tsd = 720

nstart = int(day_s*86400/dumpfreq) # integer of dumpfrec, either start from zero or 
nend = int(day_e*86400/dumpfreq) # no of time step 

f0 = 6.5e-5
rho0 = 999.8
Rmax = 25e3
vmax = 0.26

Ro=vmax/(f0*Rmax)

# Unit for the plot
conv_factor = 86400 # second in day to convert to m/day
conv_unit = 'm/day'

# Define time step 
itrs  = np.arange(nstart,nend,1)
itrsd = np.arange(day_s,day_e,1)
 
nit =itrs.size
nitd= itrsd.size

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')

nr = RC.size
nx = XC[0,:].size
ny = YC[:,0].size
icx = int(nx/2)
icy = int(ny/2)

XG = mit.rdmds('XG')
YG = mit.rdmds('YG')
dXC = mit.rdmds('DXC')
dYC = mit.rdmds('DYC')
dXG = mit.rdmds('DXG')
dYG = mit.rdmds('DYG')

xc_dom =XC[plta:pltb, plta:pltb]*1e-3
yc_dom =YC[plta:pltb, plta:pltb]*1e-3

Vort_c = np.zeros((ny,nx));
Vort_all = np.zeros((nit,ny,nx));

termu = np.zeros((ny,nx));
termv = np.zeros((ny,nx));
#
dytermu = np.zeros((ny,nx));
dxtermv = np.zeros((ny,nx));
W_pall = np.zeros((nitd,ny,nx));

W5d=np.zeros((nit,nr,ny,nx))

idepth = 65 # z indice
int_depth = 25

dir0 = './'
#./../MITgcm/mitgcm_configs/eddy_airsea/run51a_bc_a25_ua_Av1e-2'

file1 = 'diagU*'
file2 = 'diagV*'
file3 = 'diagSurf*'

# time loop
for it in itrs:
    U = mit.rdmds('U',((it)*ts))
    V = mit.rdmds('V',((it)*ts))
    W = mit.rdmds('W',((it)*ts))
    Vort = vorticity(U[0,:,:],V[0,:,:],dXC,dYC)
    W5d[int((it-itrs[0])), :,:,:]= W;
#
#Interpolate to cell center
    for z in np.arange(0,nr):
        for y in np.arange(0,ny):
            Vort_c[y,:] = np.interp(XC[y,:],XG[y,:],Vort[y,:])
        for x in np.arange(0,nx):
            Vort_c[:,x] = np.interp(YC[:,x],YG[:,x],Vort_c[:,x]) #verify if V points = YUg

    Vort_all[int((it-itrs[0])),:,:]= Vort_c;
#
Vort_av=np.mean(Vort_all, axis =0)

for it in itrsd:
#
    upress = mit.rdmds(dir0 + file1,it*tsd,rec=5)
    vpress = mit.rdmds(dir0 + file2,it*tsd,rec=5)
    int_upress = np.trapz(upress[0:int_depth,:,:],-RC[0:int_depth].squeeze(), axis=0)
    int_vpress = np.trapz(vpress[0:int_depth,:,:],-RC[0:int_depth].squeeze(), axis=0)
    termu = int_upress/(f0+Vort_av)
    termv = int_vpress/(f0+Vort_av)
#
    for i in range(1,nx-2):
        for j in range(1,ny-2):
            dytermu[(j+1),i] = (termu[j+1,i]-termu[j,i])/dYC[j+1,i];
            dxtermv[j,(i+1)] = (termv[j,i+1]-termv[j,i])/dXC[j,i+1];
#            W_p[(j+1),(i+1)] = (dxtermv[j+1,i+1]-dytermu[j+1,i+1]);
    W_p = dxtermv-dytermu
    W_pall[int((it-itrsd[0])),:,:]= W_p; #u-point u_wt

W_ta = np.mean(W5d, axis=0)
W_pa = np.mean(W_pall, axis=0)

# ==================== plot 

we_range = np.linspace(-5.0, 5.0, 101, endpoint=True)
we_ticks = np.linspace(-5.0, 5.0, 11, endpoint=True)
fig = plt.figure(figsize=(10,5))
fig.suptitle('Averaged Vertical Velocity, day %d to %d' % (day_s, day_e),fontsize=12)
#
ax1 = plt.subplot(121)
plt.contourf(xc_dom,yc_dom,W_ta[int_depth,plta:pltb,plta:pltb]*conv_factor,we_range,cmap=cm.seismic)#
cb=plt.colorbar(ticks=we_ticks, format='%1.1f')
cb.set_label(r'$W \, (m/day)$', labelpad=-40, y=1.1, rotation=0)
rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
ax1.add_artist(rmax1)
rmax12 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax1.add_artist(rmax12)
ax1.set_aspect(1)
plt.text(posext,posexty,'$W_{extrema}=$ [%1.1f, %1.1f] $m/day$' % (np.min(W_ta[int_depth,:,:])*conv_factor,np.max(W_ta[int_depth,:,:])*conv_factor), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'Time averaged, $\overline{W}$ at $z=%d$' % (RC[int_depth]),fontsize=11)

ax1 = plt.subplot(122)
plt.contourf(xc_dom,yc_dom,W_pa[plta:pltb,plta:pltb]*conv_factor,we_range,cmap=cm.seismic)#
cb=plt.colorbar(ticks=we_ticks, format='%1.1f')
cb.set_label(r'$W \, (m/day)$', labelpad=-40, y=1.1, rotation=0)
rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
ax1.add_artist(rmax1)
rmax12 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax1.add_artist(rmax12)
ax1.set_aspect(1)
plt.text(posext,posexty,'$W_{extrema}=$ [%1.1f, %1.1f] $m/day$' % (np.min(W_pa[:,:])*conv_factor,np.max(W_pa[:,:])*conv_factor), fontsize=10)
plt.xlabel("x (km)")
#plt.ylabel("y (km)")
plt.title('W from pressure gradient', fontsize=11)
plt.tight_layout (pad = 1)
plt.savefig('./figures/Wpress_grad_day%d-%d.png' % (day_s, day_e))

