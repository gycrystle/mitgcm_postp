"""
Compute 
W_ekman : model produced Ekman pumping from ek transport divergence
W_stern : Theoretical nonlinear Ekman pumping, Stern, 1965
W_wt    : Theoretical nonlinear Ekman pumping w/ curvature, Wenegrat & Thomas, 2017

"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from scipy import interpolate
from matplotlib import cm
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

def radial_prof(data, r):
    uniq = np.unique(r)
    prof = np.array([ np.mean(data[ r==un ]) for un in uniq ])
    return uniq, prof

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
plta = 30
pltb = 150
idepth = 99 #depth to substract geostrophic flow 
posext = -100 # location of text 

timestep = 180 #timestep input in second
dumpfreq = 10800 # in second
ts = dumpfreq/timestep  # time step = ts*dt (in second); = 7200 = dumpfreq

day_s = 10
day_e = 15

nstart = int(day_s*86400/dumpfreq) # integer of dumpfrec, either start from zero or 
nend = int(day_e*86400/dumpfreq) # no of time step 
#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(nstart,nend,1)
nit=itrs.size

itpd = int(86400/dumpfreq)

iters1 = mit.mds.scanforfiles('diagTAUX')

iters2=np.intersect1d(itrs*ts, iters1)
nitt = iters2.size
itrst=np.arange(0,nitt, 1)
tst=iters2[1]-iters2[0]

itrsd=np.arange(day_s+1,day_e+1,1)
nitd=itrsd.size
tsd = 480


f0 = 6.5e-5
rho0 = 999.8
Rmax = 25e3
vmax = 0.26
Ro=vmax/(f0*Rmax)

conv_factor = 86400 # second in day to convert to m/day
conv_unit = 'm/day'

# Load data
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

dX_Ug = np.tile(np.diff(XC, axis=1),[nr,1,1])
dY_Vg = np.tile(np.diff(YC, axis=0),[nr,1,1])

xc_dom =XC[plta:pltb, plta:pltb]*1e-3-XC[icy,icx]*1e-3
yc_dom =YC[plta:pltb, plta:pltb]*1e-3-YC[icy,icx]*1e-3

xc_domo =XC[plta:pltb, plta:pltb]*1e-3
yc_domo =YC[plta:pltb, plta:pltb]*1e-3


YUg = YC+0.5*dYC # grid where Ug calculated
XUg = XC+0.5*dXC # grid where Vg calculated

Ug_all=np.zeros((nit,nr,ny,nx))
Vg_all=np.zeros((nit,nr,ny,nx))

Utheta_all=np.zeros((nit,nr,ny,nx))
Ut_all=np.zeros((nit,nr,ny,nx))
Vt_all=np.zeros((nit,nr,ny,nx))

u_wtall=np.zeros((nit,nr,ny,nx))
v_wtall=np.zeros((nit,nr,ny,nx))
u_ekall=np.zeros((nit,nr,ny,nx))
v_ekall=np.zeros((nit,nr,ny,nx))

Ug_c=np.zeros((nr,ny,nx))
Vg_c=np.zeros((nr,ny,nx))
Ug300_c=np.zeros((ny,nx))
Vg300_c=np.zeros((ny,nx))

Ut_c=np.zeros((nr,ny,nx))
Vt_c=np.zeros((nr,ny,nx))
Wt_c=np.zeros((nr,ny,nx))
Vort_g =np.zeros((ny,nx))
u_ek = np.zeros((nr,ny,nx))
v_ek = np.zeros((nr,ny,nx))

# To compute W_stern
dxU = np.zeros((ny,nx));
dxV = np.zeros((ny,nx));
dxU_c = np.zeros((ny,nx));
dxV_c = np.zeros((ny,nx));


zeta = np.zeros((ny,nx));
zeta_intx = np.zeros((ny,nx));
zeta_inty = np.zeros((ny,nx));

termx = np.zeros((ny,nx));
termy = np.zeros((ny,nx));
#
dytermx = np.zeros((ny,nx));
dxtermy = np.zeros((ny,nx));
W_stern = np.zeros((ny,nx));

Vort_all = np.zeros((nit,ny,nx))
Vort_da = np.zeros((nitt,ny,nx));
U_da = np.zeros((nitt,nr,ny,nx));
V_da = np.zeros((nitt,nr,ny,nx));
Utheta_da = np.zeros((nitt,nr,ny,nx));
W5d = np.zeros((nit,nr,ny,nx))

#       transform to polar coord
thetaC = np.arctan2(YC-YC[icy,icx], XC-XC[icy,icx])
radC = np.hypot(XC-XC[icy,icx],YC-YC[icy,icx])

# time loop
for it in itrs:
    U = mit.rdmds('U',((it+1)*ts))
    V = mit.rdmds('V',((it+1)*ts))
    Vort = vorticity(U[99,:,:],V[99,:,:],dXG,dYG)

    W = mit.rdmds('W',((it+1)*ts))
#    W5d[int((it-itrs[0])), :,:,:]= W;
# === To compute U ageostrophic from U bottom
    u_eki = U - U[idepth,:,:]
    v_eki = V - V[idepth,:,:]
# ========= trial
#    Vort_use = Vort_da[it,:,:]#Vort300 # or Vort_da[it-itrsd[0],:,:]
    Vort_f0 = Vort[plta:pltb,plta:pltb]/f0 #or Vort_da[it-itrsd[0],plta:pltb,plta:pltb]/f0
    zetamin = np.min(Vort_f0)
    zetamax = np.max(Vort_f0)
    if (abs(zetamin)>zetamax):
        vort_contour = [zetamin*2/3]
    else:
        vort_contour = [zetamax*2/3]
#
#       define eddy center
    plt.figure()
    CS = plt.contour(xc_domo,yc_domo,Vort_f0, vort_contour)
    c =  center_of_mass(CS.allsegs[0][0])
    plt.close()

    XGi = XG-c[0]*1e3+XG[icy,icx]
    YGi = YG-c[1]*1e3+YG[icy,icx]
    # Interpolate to cell center
#    fvort = interpolate.interp2d(XGi,YGi,Vort, kind='cubic')
#    Vort_g[:,:] = fvort(XG[icy,:],YG[:,icx])

    for z in np.arange(0,nr):
#        fU = interpolate.interp2d(XGi,YGi,U[z,:,:], kind='cubic')
#        Ut_c[z,:,:] = fU(XC[icy,:],YC[:,icx])
#        fV = interpolate.interp2d(XGi,YGi,V[z,:,:], kind='cubic')
#        Vt_c[z,:,:] = fV(XC[icy,:],YC[:,icx])
#        fuek = interpolate.Rbf(XGi[icy,:],YGi[:,icx],u_eki[z,:,:], function='multiquadric')
#        u_ek[z,:,:] = fuek(XG,YC)
#        fvek = interpolate.Rbf(XGi[icy,:],YGi[:,icx],v_eki[z,:,:], function='multiquadric')
#        v_ek[z,:,:] = fvek(XC,YG) #np.interp(XG[y,:],XGi[y,:],v_eki[z,y,:])
#        fW = interpolate.Rbf(XGi[icy,:],YGi[:,icx],W[z,:,:], function='multiquadric')
#        Wt_c[z,:,:] = fW(XC,YC)

        for y in np.arange(0,ny):
            Ut_c[z,y,:] = np.interp(XC[y,:],XGi[y,:],U[z,y,:])
            Vt_c[z,y,:] = np.interp(XC[y,:],XGi[y,:],V[z,y,:])
            u_ek[z,y,:] = np.interp(XG[y,:],XGi[y,:],u_eki[z,y,:])
            v_ek[z,y,:] = np.interp(XC[y,:],XGi[y,:],v_eki[z,y,:])
            Vort_g[y,:] = np.interp(XG[y,:],XGi[y,:],Vort[y,:])
            Wt_c[z,y,:] = np.interp(XC[y,:],XGi[y,:],W[z,y,:])
        for x in np.arange(0,nx):
            Vt_c[z,:,x] = np.interp(YC[:,x],YGi[:,x],Vt_c[z,:,x]) #verify if V points = YUg
            Ut_c[z,:,x] = np.interp(YC[:,x],YGi[:,x],Ut_c[z,:,x])
            u_ek[z,:,x] = np.interp(YC[:,x],YGi[:,x],u_ek[z,:,x]) #verify if V points = YUg
            v_ek[z,:,x] = np.interp(YG[:,x],YGi[:,x],v_ek[z,:,x])
            Vort_g[:,x] = np.interp(YG[:,x],YGi[:,x],Vort_g[:,x])
            Wt_c[z,:,x] = np.interp(YC[:,x],YGi[:,x],Wt_c[z,:,x])

    Ut_all[int((it-itrs[0])), :,:,:] = Ut_c
    Vt_all[int((it-itrs[0])), :,:,:] = Vt_c
    Vort_all[int((it-itrs[0])), :,:] = Vort_g
    W5d[int((it-itrs[0])), :,:,:]= Wt_c;
    u_ekall[(it-itrs[0]),:,:,:]= u_ek; #u-point u_wt
    v_ekall[(it-itrs[0]),:,:,:]= v_ek; #v-point v_wti

#       transform to polar coord
#    thetaC = np.arctan2(YC-c[1]*1e3, XC-c[0]*1e3)
#    radC = np.hypot(XC-c[0]*1e3,YC-c[1]*1e3)

#       Compute Tangential and Angular velocity
    U_theta  = Vt_c*np.cos(thetaC) - Ut_c*np.sin(thetaC)
    Utheta_all[int((it-itrs[0])), :,:,:] = U_theta
#
# ============== trial end
#
#
#    if (((it+1) % itpd)==0):# and (it-itrs[0] != 0):
    if ((it+1)*ts) in iters1:# and (it-itrs[0] != 0):
        itnitt = int((it+1-itrs[0])/(tst/ts))-1
        Vort_da[itnitt,:,:]=np.nanmean(Vort_all[(it-itrs[0]-3):(it+1-itrs[0]),:,:], axis=0)
        U_da[itnitt,:,:,:]=np.nanmean(Ut_all[(it-itrs[0]-3):(it+1-itrs[0]),:,:], axis=0)
        V_da[itnitt,:,:,:]=np.nanmean(Vt_all[(it-itrs[0]-3):(it+1-itrs[0]),:,:], axis=0)
        Utheta_da[itnitt,:,:,:]=np.nanmean(Utheta_all[(it-itrs[0]-3):(it+1-itrs[0]),:,:], axis=0)

iters2=np.intersect1d(itrs*ts, iters1)
nitt = iters2.size
itrst=np.arange(0,nitt, 1)
tst=iters2[1]-iters2[0]

itrsd=np.arange(day_s+1,day_e+1,1)
nitd=itrsd.size
tsd = 480

W_sternall = np.zeros((nitt,ny,nx));
W_b = 0.0*W_sternall
W_a = 0.0*W_sternall
Mx_WTall = np.zeros((nitt,ny,nx));
My_WTall = np.zeros((nitt,ny,nx));
Vort_cndall = np.zeros((nitt,ny,nx));
Omega_all = np.zeros((nitt,ny,nx));

Vort_c= np.zeros((ny,nx));

for it in itrst:
# ----------------------- #1 Stern -----------------------------------
#    if (it % itpd)==0:
    U = mit.rdmds('U',(iters2[it]))
    V = mit.rdmds('V',(iters2[it]))
#    Vort = vorticity(U[99,:,:],V[99,:,:],dXC,dYC)

    taux = mit.rdmds('diagTAUX',(iters2[it]))
    tauy = mit.rdmds('diagTAUY',(iters2[it]))
#    Vort_da=np.nanmean(Vort_all[(it-itrsd[0]-8):(it-itrs[0]),:,:], axis=0)
#
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            zeta_intx[j,i] =(Vort_da[it,j,i]+Vort_da[it,j+1,i])/2 #vorticity @ u point
            zeta_inty[j,i] =(Vort_da[it,j,i]+Vort_da[it,j,i+1])/2
#
            termx[j,i] =taux[j,i]/(f0+zeta_intx[j,i])
            termy[j,i] =tauy[j,i]/(f0+zeta_inty[j,i])
#
    for i in range(1,nx-2):
        for j in range(1,ny-2):
            dytermx[(j+1),i] = (termx[j+1,i]-termx[j,i])/dYC[j+1,i];
            dxtermy[j,(i+1)] = (termy[j,i+1]-termy[j,i])/dXC[j,i+1];
            W_stern[(j+1),(i+1)] = (dxtermy[j+1,i+1]-dytermx[j+1,i+1])/rho0;
#
#    W_sternall[int((it-itrs[0])/itpd),:,:]= W_stern; #u-point u_wt
    W_sternall[it,:,:]= W_stern; #u-point u_wt
# ------------------------------------------------------------------- 
# --------------------  Wenegrat & Thomas ---------------------------
#    Vort300 = vorticity(U[99,:,:],V[99,:,:],dXC,dYC)
    Vort_use = Vort_da[it,:,:]#Vort300 # or Vort_da[it-itrsd[0],:,:]
    """
    Vort_f0 = Vort_use[plta:pltb,plta:pltb]/f0 #or Vort_da[it-itrsd[0],plta:pltb,plta:pltb]/f0
    zetamin = np.min(Vort_f0)
    zetamax = np.max(Vort_f0)
    if (abs(zetamin)>zetamax):
        vort_contour = [zetamin*2/3]
    else:
        vort_contour = [zetamax*2/3]
#
#       define eddy center
    plt.figure()
    CS = plt.contour(xc_dom,yc_dom,Vort_f0, vort_contour)
    c =  center_of_mass(CS.allsegs[0][0])
    plt.close()

    # Interpolate to cell center, loosing first grids(bottom and left)
    for z in np.arange(0,nr):
        for y in np.arange(0,ny):
            Ut_c[z,y,:] = np.interp(XC[y,:],XG[y,:],U_da[it,z,y,:])
        for x in np.arange(0,nx):
            Vt_c[z,:,x] = np.interp(YC[:,x],YG[:,x],V_da[it,z,:,x]) #verify if V points = YUg

#       transform to polar coord
    thetaC = np.arctan2(YC-c[1]*1e3, XC-c[0]*1e3)
    radC = np.hypot(XC-c[0]*1e3,YC-c[1]*1e3)
    """
#       Compute Tangential and Angular velocity
    U_theta = Utheta_da[it] #Vt_c*np.cos(thetaC) - Ut_c*np.sin(thetaC)
# Compute azimuthal average
    [r_av, Utheta_av] = radial_prof(U_theta[99,:,:], radC)
    Omega_av = Utheta_av/r_av
    drOmega_av = Omega_av*0.0
    for j in range(1,r_av.size):
        drOmega_av[j] = (Omega_av[j]-Omega_av[j-1])/(r_av[j]-r_av[j-1])#/np.diff(r_av)
#projecting back to grid
    U_theta_av=U_theta[99,:,:]*0.0
    dr_Omega_av = U_theta[99,:,:]*0.0
    for ir in range(0,r_av.size):
        indicer= np.where(radC==r_av[ir])
#        for indi in range(0,indicer[0].size): 
        U_theta_av[indicer[0][:],indicer[1][:]] = Utheta_av[ir]
        dr_Omega_av[indicer[0][:],indicer[1][:]] = drOmega_av[ir]
#    Omega = U_theta[99,:,:]/radC/(0.26/(25*1e3))
    Omega = U_theta_av/radC/(0.26/(25*1e3))

    for y in np.arange(0,ny):
        Vort_c[y,:] = np.interp(XC[y,:],XG[y,:],Vort_use[y,:])
    for x in np.arange(0,nx):
        Vort_c[:,x] = np.interp(YC[:,x],YG[:,x],Vort_c[:,x]) #verify if V points = YUg
    Vort_cnd = Vort_c/(0.26/(25*1e3))

#   Compute Ekman transport Wenegrat & Thomas
    dim_factor = np.max(taux)/(rho0*f0)
    Mx_WT = dim_factor*Ro*(Vort_cnd-2*Omega)*np.sin(thetaC)*np.cos(thetaC)/((1+Ro*2*Omega)*(1+Ro*Vort_cnd)-Ro**2*Omega**2)
    My_WT = -dim_factor*(1+Ro*Omega+Ro*2*Omega*np.sin(thetaC)**2+Ro*Vort_cnd*np.cos(thetaC)**2)/((1+Ro*2*Omega)*(1+Ro*Vort_cnd)-Ro**2*Omega**2)
#     Small Ro
#    Mx_WT = (np.max(taux)/(rho0*f0))*Ro*(Vort_cnd-2*Omega)*np.sin(thetaC)*np.cos(thetaC)
#    My_WT = -(np.max(taux)/(rho0*f0))*(1-Ro*((Vort_cnd-Omega)*np.sin(thetaC)**2 + Omega*np.cos(thetaC)**2))
#
    Mx_WTall[it,:,:]=Mx_WT # @cell center
    My_WTall[it,:,:]=My_WT # @cell center
# -------------------------------------------------------------------
#       Compute correction term by bruno
    """
    for i in range(1,nx):
        for j in range(1,ny):
            dxV[j,i] = (V_da[it,99,j,i]-V_da[it,99,j,i-1])/dXC[j,i];
            dxU[j,i] = (U_da[it,99,j,i]-U_da[it,99,j,i-1])/dXC[j,i];

#   Interpolate dxV from cell corner to cell center
#    dxV_c = 0.0*dxV
    for y in np.arange(0,ny):
        dxV_c[y,:] = np.interp(XC[y,:],XUg[y,:],dxV[y,:])
        dxU_c[y,:] = np.interp(XC[y,:],XUg[y,:],dxU[y,:])
#    for x in np.arange(0,nx):
#        dxV_c[:,x] = np.interp(YC[:,x],YG[:,x],dxV_c[:,x])
#
    drUth=-dxU_c*np.sin(thetaC)*np.cos(thetaC)+dxV_c*(np.cos(thetaC)**2)
    drOmega=(drUth-U_theta_av/radC)/radC
    omega_0 = (2/radC-1)*U_theta[0,:,:] + drUth #vorticity calculated in cylindrical coord
    """
#
    ABC = taux*np.sin(thetaC)*dr_Omega_av/rho0
    B = ABC/((f0+Vort_c)**2) #correction term bruno
    AS = ABC/((f0+Vort_c)*(f0+2*U_theta_av/radC)) #correction term alex
#    B = taux*np.sin(thetaC)*(drUth-U_theta_av/radC)/radC/(f0+Vort_c)**2/rho0 #correction term bruno
#    B = taux*np.sin(thetaC)*(drUth-U_theta[99,:,:]/radC)/radC/(f0+Vort_c)**2/rho0 #correction term bruno
##    B = taux*np.sin(thetaC)*(drUth-U_theta[0,:,:]/radC)/radC/(f0+omega_0)**2/rho0 #correction term bruno
#    AS = taux*np.sin(thetaC)*(drUth-U_theta_av/radC)/radC/((f0+Vort_c)*(f0+2*omega_0)*rho0) #correction term bruno
    W_b[it,:,:]= W_stern+B;
    W_a[it,:,:]= W_stern+AS;


#======= time loop end
#
# Compute ekman transport
Mx_ek = np.trapz(u_ekall[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1)
My_ek = np.trapz(v_ekall[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1)

# Compute divergence of Ekman transport
# Wenegrat & Thomas
dxMx_WT = np.diff(Mx_WTall, axis=2)/np.tile(dXC[:,1:],[nitt,1,1]) # V-point
dyMy_WT = np.diff(My_WTall, axis=1)/np.tile(dYC[1:,:],[nitt,1,1]) # V-point
# Model produced
dxMx_ek = np.diff(Mx_ek, axis=2)/np.tile(dXC[:,1:],[nit,1,1]) # cell center
dyMy_ek = np.diff(My_ek, axis=1)/np.tile(dYC[1:,:],[nit,1,1]) # cell center

#interpolate to cell center
nitd = Mx_WTall.shape[0]
dxMx_WTi=np.zeros((nitt,ny,nx))
dyMy_WTi=np.zeros((nitt,ny,nx))
for it in np.arange(0,nitt):
    for y in np.arange(0,ny):
        dxMx_WTi[it,y,:] = np.interp(XC[y,:],XG[y,1:],dxMx_WT[it,y,:])
    for x in np.arange(0,nx):
        dyMy_WTi[it,:,x] = np.interp(YC[:,x],YG[1:,x],dyMy_WT[it,:,x])

W_ek = dxMx_ek[:,1:,:] + dyMy_ek[:,:,1:] # no need to interpolate because already at cell center
W_WTi = dxMx_WTi + dyMy_WTi

Mx_WTm = np.mean(Mx_WTall, axis=0)
My_WTm = np.mean(My_WTall, axis=0)

Mx_ekman = np.mean(Mx_ek, axis=0)#/(tau0*Roe/(rho0*f0))
My_ekman = np.mean(My_ek, axis=0)#/(tau0/(rho0*f0))

# === Time averaging all Ekman pumping ===

# W time averaged
W_ta = np.mean(W5d, axis=0)

# Model produced Ekman transport div
W_ekman = np.mean(W_ek, axis=0) #/(tau0*Roe/(rho0*f0*Rmax))

# Non-linear Ekman pumping (Stern)
Wm_stern=np.mean(W_sternall, axis=0)

# Ekman pumping Wenegrat & Thomas
W_WT = np.mean(W_WTi, axis=0)#/(tau0*Roe/(rho0*f0*Rmax))

We_bd = np.mean(W_b, axis=0)
We_as = np.mean(W_a, axis=0)
#
###  ============ PLOT ============
# W_ekman bottom 
#
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
#
time = (itrs)*dumpfreq/3600
maxdepthi = np.where(W_ta==np.max(W_ta))
maxdepth = maxdepthi[0][0]
maxdepth = 17
#
def subplot_plan(W,title):
    plt.contourf(xc_dom,yc_dom,W[maxdepth,plta:pltb,plta:pltb]*conv_factor,we_range,cmap=cm.seismic)#
    cb=plt.colorbar(ticks=we_ticks, format='%1.1f')
    cb.set_label(r'$W \, (m/day)$', labelpad=-40, y=1.1, rotation=0)
    rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
    ax1.add_artist(rmax1)
    rmax12 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
    ax1.add_artist(rmax12)
    #CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
    #plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
    ax1.set_aspect(1)
    plt.text(posext,posext,'$W_{extrema}=$ [%1.1f, %1.1f] $m/day$' % (np.min(W[maxdepth,:,:])*conv_factor,np.max(W[maxdepth,:,:])*conv_factor), fontsize=10)
    #plt.xlabel("x (km)")
    plt.ylabel("y (km)")
    plt.title(title,fontsize=11)
#
## ===================  Plot Ekman pumping ===============================
we_range = np.linspace(-6.5, 6.5, 101, endpoint=True)
we_ticks = np.linspace(-6.5, 6.5, 11, endpoint=True)
fig = plt.figure(figsize=(15,8))
fig.suptitle('Ekman pumping, day %d to %d' % (day_s, day_e),fontsize=12)
#
ax1 = plt.subplot(231)
title1="r'$\overline{W_{e}}=\overline{\partial_x \int \, U_e + \partial_y \int \, V_e}$'"
subplot_plan(W_ta,title1)
"""
plt.contourf(xc_dom,yc_dom,W_ta[maxdepth,plta:pltb,plta:pltb]*conv_factor,we_range,cmap=cm.seismic)#
cb=plt.colorbar(ticks=we_ticks, format='%1.1f')
cb.set_label(r'$W \, (m/day)$', labelpad=-40, y=1.1, rotation=0)
rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
ax1.add_artist(rmax1)
rmax12 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax1.add_artist(rmax12)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
ax1.set_aspect(1)
plt.text(posext,posext,'$W_{extrema}=$ [%1.1f, %1.1f] $m/day$' % (np.min(W_ta[maxdepth,:,:])*conv_factor,np.max(W_ta[maxdepth,:,:])*conv_factor), fontsize=10)
#plt.text(posext,posext,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(W_Bm)*1e3,np.max(W_Bm)*1e3), fontsize=10)
#plt.xlabel("x (km)")
plt.ylabel("y (km)")
#plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
plt.title(r'Time averaged, $\overline{W}$ at $z=%d$' % (RC[maxdepth]),fontsize=11)
"""
# ===============================================
ax2 = plt.subplot(232, sharex=ax1)
title2="r'$\overline{W_{e}}=\overline{\partial_x \int \, U_e + \partial_y \int \, V_e}$'"
subplot_plan(W_ekman,title2)
"""
plt.contourf(xc_dom,yc_dom,W_ekman[plta:pltb,plta:pltb]*conv_factor,we_range,cmap=cm.seismic)#
cb=plt.colorbar(ticks=we_ticks, format='%1.1f')
cb.set_label(r'$W \, (m/day)$', labelpad=-40, y=1.1, rotation=0)
rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
ax2.add_artist(rmax1)
rmax12 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax2.add_artist(rmax12)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
ax2.set_aspect(1)
plt.text(posext,posext,'$W_{extrema}=$ [%1.1f, %1.1f] $m/day$' % (np.min(W_ekman)*conv_factor,np.max(W_ekman)*conv_factor), fontsize=10)
#plt.xlabel("x (km)")
#plt.ylabel("y (km)")
#plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
plt.title(r'$\overline{W_{e}}=\overline{\partial_x \int \, U_e + \partial_y \int \, V_e}$',fontsize=11)
"""
# ===============================================
#
ax3 = plt.subplot(233, sharex=ax1)
title3="r'$\overline{W_{stern}}=\overline{\frac{1}{\rho_0} \nabla \times \left[\frac{\tau}{f+\zeta}\right]}$'"
subplot_plan(Wm_stern,title3)
ax3.set_aspect(1)
"""
im=plt.contourf(xc_dom,yc_dom,Wm_stern[plta:pltb,plta:pltb]*conv_factor,we_range,cmap=cm.seismic)#
#divider = make_axes_locatable(plt.gca())
#cax = divider.append_axes("right", "5%", pad="3%")
#cb=plt.colorbar(im,ticks=we_ticks, format='%1.2f', cax=cax)

cb=plt.colorbar(ticks=we_ticks, format='%1.1f')
cb.set_label(r'$W \, (m/day)$', labelpad=-40, y=1.1, rotation=0)
rmax2 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
ax3.add_artist(rmax2)
rmax22 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax3.add_artist(rmax22)
#CS1 = plt.contour(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3, levels, colors='0.6')
#plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)
plt.text(posext,posext,'$W_{extrema}=$ [%1.1f, %1.1f] $m/day$' % (np.min(Wm_stern)*conv_factor,np.max(Wm_stern)*conv_factor), fontsize=10)
#plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
plt.title(r'$\overline{W_{stern}}=\overline{\frac{1}{\rho_0} \nabla \times \left[\frac{\tau}{f+\zeta}\right]}$', fontsize=14)
"""
# ================================================
#
ax4 = plt.subplot(234, sharex=ax1)
title4 = r'$\overline{W_{wt}}$, Wenegrat&Thomas'
subplot_plan(W_WT,title4)
"""
plt.contourf(xc_dom,yc_dom,W_WT[plta:pltb,plta:pltb]*conv_factor,we_range,cmap=cm.seismic)#
#plt.contourf(xc_dom,yc_dom,W_Bm[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(ticks=we_ticks, format='%1.1f')
cb.set_label(r'$W \, (m/day)$', labelpad=-40, y=1.1, rotation=0)
rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
ax4.add_artist(rmax1)
rmax12 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax4.add_artist(rmax12)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
ax4.set_aspect(1)
plt.text(posext,posext,'$W_{extrema}=$ [%1.1f, %1.1f] $m/day$' % (np.nanmin(W_WT)*conv_factor,np.nanmax(W_WT)*conv_factor), fontsize=10)
#plt.text(posext,posext,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(W_Bm)*1e3,np.max(W_Bm)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
plt.title(r'$\overline{W_{wt}}$, Wenegrat&Thomas',fontsize=10)
"""
# ===============================================
# ================================================
#
ax5 = plt.subplot(235, sharex=ax1)
title5 = r'$\overline{W_{b}}=W_{stern}+ \frac{\partial_r \Omega}{\rho (f+\zeta )^2}$, Bruno calculation'
subplot_plan(We_bd,title5)
"""
plt.contourf(xc_dom,yc_dom,We_bd[plta:pltb,plta:pltb]*conv_factor,we_range,cmap=cm.seismic)#
#plt.contourf(xc_dom,yc_dom,W_Bm[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(ticks=we_ticks, format='%1.1f')
cb.set_label(r'$W \, (m/day)$', labelpad=-40, y=1.1, rotation=0)
rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
ax5.add_artist(rmax1)
rmax12 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax5.add_artist(rmax12)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
ax5.set_aspect(1)
plt.text(posext,posext,'$W_{extrema}=$ [%1.1f, %1.1f] $m/day$' % (np.nanmin(We_bd)*conv_factor,np.nanmax(We_bd)*conv_factor), fontsize=10)
#plt.text(posext,posext,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(W_Bm)*1e3,np.max(W_Bm)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)

plt.title(r'$\overline{W_{b}}=W_{stern}+ \frac{\partial_r \Omega}{\rho (f+\zeta )^2}$, Bruno calculation',fontsize=10)
"""
# ===============================================

#
ax6 = plt.subplot(236, sharex=ax1)
title6 = r'$\overline{W_{a}}=W_{stern}+ \frac{\partial_r \Omega}{\rho (f+\zeta ) (f+2\Omega )}$, Alex correction'
subplot_plan(We_as, title6)
# ================================================
#
"""
plt.contourf(xc_dom,yc_dom,We_as[plta:pltb,plta:pltb]*conv_factor,we_range,cmap=cm.seismic)#
#plt.contourf(xc_dom,yc_dom,W_Bm[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(ticks=we_ticks, format='%1.1f')
cb.set_label(r'$W \, (m/day)$', labelpad=-40, y=1.1, rotation=0)
rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
ax5.add_artist(rmax1)
rmax12 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax5.add_artist(rmax12)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
ax5.set_aspect(1)
plt.text(posext,posext,'$W_{extrema}=$ [%1.1f, %1.1f] $m/day$' % (np.nanmin(We_as)*conv_factor,np.nanmax(We_as)*conv_factor), fontsize=10)
#plt.text(posext,posext,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(W_Bm)*1e3,np.max(W_Bm)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
plt.title(r'$\overline{W_{a}}=W_{stern}+ \frac{\partial_r \Omega}{\rho (f+\zeta ) (f+2\Omega )}$, Alex correction',fontsize=10)
"""
# ===============================================
plt.savefig('./figures/W_ek_mday_day%d-%d.png' % (day_s, day_e))

fig = plt.figure(figsize=(5,4))
#ax5.set_aspect(1)
plt.plot(yc_dom[:,icx],W_ta[maxdepth,plta:pltb,icx]*conv_factor,label=r'$\overline{W}$')# label='Time averaged')
plt.plot(yc_dom[:,icx],W_ekman[plta:pltb,icx]*conv_factor, label=r'$\overline{W_{e}}$') #r'$\overline{\partial_x \int \, U_e + \partial_y \int \, V_e}$')
plt.plot(yc_dom[:,icx],Wm_stern[plta:pltb,icx]*conv_factor, label=r'$\overline{W_{stern}}$')
plt.plot(yc_dom[:,icx],W_WT[plta:pltb,icx]*conv_factor, label=r'$\overline{W_{wt}}$')
plt.plot(yc_dom[:,icx],We_bd[plta:pltb,icx]*conv_factor, label=r'$\overline{W_b}$')
plt.plot(yc_dom[:,icx],We_as[plta:pltb,icx]*conv_factor, label=r'$\overline{W_a}$')
plt.ylim(-8,8)
plt.legend(loc='lower right', fancybox=True, framealpha=0.5)
plt.grid()
plt.xlabel("y (km)")
plt.ylabel("Vertical velocity (m/day)")
#plt.title("Comparison",fontsize=10)
#
#plt.tight_layout(pad=1)
plt.savefig('./figures/W_ek_compare_day%d-%d.png' % (day_s, day_e))

## ================= PLOT Ekman Transport ====================

#we_range = np.linspace(-0.02, 0.02, 101, endpoint=True)
fig = plt.figure(figsize=(5,8))
ax1 = plt.subplot(211)
plt.contourf(xc_dom,yc_dom,Mx_ekman[plta:pltb, plta:pltb],101,norm=MidpointNormalize(midpoint=Mx_ekman[plta+1,plta+1]),cmap=cm.seismic)
cb = plt.colorbar(format='%1.3f')
cb.set_label(r'$M_e \, (m^2/s)$', labelpad=-40, y=1.1, rotation=0)
rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
rmax11 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax1.add_artist(rmax1)
ax1.add_artist(rmax11)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
ax1.set_aspect(1)
#plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(Wm_wt)*1e3,np.max(Wm_wt)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
plt.title(r'$\overline{M_e^x}=\overline{\int \, U_e}$ day %d-%d' % (day_s, day_e),fontsize=11)
#
ax2 = plt.subplot(212, sharex=ax1)
ax2.set_aspect(1)
plt.contourf(xc_dom,yc_dom,My_ekman[plta:pltb, plta:pltb],101,norm=MidpointNormalize(midpoint=My_ekman[pltb,pltb]),cmap=cm.seismic)
cb=plt.colorbar(format='%1.3f')
cb.set_label(r'$M_e \, (m^2/s)$', labelpad=-40, y=1.1, rotation=0)
rmax2 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
rmax12 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax2.add_artist(rmax2)
ax2.add_artist(rmax12)
#CS1 = plt.contour(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3, levels, colors='0.6')
#plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)
#plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(Wm_stern)*1e3,np.max(Wm_stern)*1e3), fontsize=10)

plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
plt.title(r'$\overline{M_e^y}=\overline{\int \, V_e}$ day %d-%d' % (day_s, day_e),fontsize=11)
#plt.text(520,-320,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
plt.tight_layout(pad = 1)
plt.savefig('./figures/Me_day%d-%d.png' % (day_s, day_e))
#    plt.close()
#
#==================== Ekman transport WT
fig = plt.figure(figsize=(5,8))
ax1 = plt.subplot(211)
plt.contourf(xc_dom,yc_dom,Mx_WTm[plta:pltb, plta:pltb],101,norm=MidpointNormalize(midpoint=Mx_WTm[plta+1,plta+1]),cmap=cm.seismic)
cb = plt.colorbar(format='%1.3f')
cb.set_label(r'$M_e \, (m^2/s)$', labelpad=-40, y=1.1, rotation=0)
rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
rmax11 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax1.add_artist(rmax1)
ax1.add_artist(rmax11)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
ax1.set_aspect(1)
#plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(Wm_wt)*1e3,np.max(Wm_wt)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
plt.title(r'$\overline{M_e^x}=\overline{\int \, U_e}$ day %d-%d' % (day_s, day_e),fontsize=11)
#
ax2 = plt.subplot(212, sharex=ax1)
ax2.set_aspect(1)
plt.contourf(xc_dom,yc_dom,My_WTm[plta:pltb, plta:pltb],101,norm=MidpointNormalize(midpoint=My_WTm[pltb,pltb]),cmap=cm.seismic)
cb=plt.colorbar(format='%1.3f')
cb.set_label(r'$M_e \, (m^2/s)$', labelpad=-40, y=1.1, rotation=0)
rmax2 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
rmax12 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax2.add_artist(rmax2)
ax2.add_artist(rmax12)
#CS1 = plt.contour(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3, levels, colors='0.6')
#plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)
#plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(Wm_stern)*1e3,np.max(Wm_stern)*1e3), fontsize=10)

plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
plt.title(r'$\overline{M_e^y}=\overline{\int \, V_e}$ day %d-%d' % (day_s, day_e),fontsize=11)
#plt.text(520,-320,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
plt.tight_layout (pad = 1)
plt.savefig('./figures/Me_WT_day%d-%d.png' % (day_s, day_e))

#
#####################################################################
