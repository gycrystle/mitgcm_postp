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
    

# select plot domain
plta = 30
pltb = 150
idepth = 99 #depth to substract geostrophic flow 
posext = 500 # location of text 

timestep = 180 #timestep input in second
dumpfreq = 10800 # in second

day_s = 10
day_e = 15

ts = dumpfreq/timestep  # time step = ts*dt (in second); = 7200 = dumpfreq

nstart = int(day_s*86400/dumpfreq) # integer of dumpfrec, either start from zero or 
nend = int(day_e*86400/dumpfreq) # no of time step 
itpd = int(86400/dumpfreq)

f0 = 6.5e-5
rho0 = 999.8
Rmax = 25e3
vmax = 0.26
Ro=vmax/(f0*Rmax)
conv_factor = 86400 # second in day to convert to m/day
conv_unit = 'm/day'
#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(nstart,nend,1)
nit=itrs.size

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')

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

nr = RC.size
nx = XC[0,:].size
ny = YC[:,0].size
icx = int(nx/2)
icy = int(ny/2)

#transform to polar coord
#thetaC = np.arctan2(YC-YC[icy, icx], XC-XC[icy, icx])
#radC = np.hypot(XC-XC[icy, icx],YC-YC[icy, icx])

XG = mit.rdmds('XG')
YG = mit.rdmds('YG')
dXC = mit.rdmds('DXC')
dYC = mit.rdmds('DYC')
dXG = mit.rdmds('DXG')
dYG = mit.rdmds('DYG')

dX_Ug = np.tile(np.diff(XC, axis=1),[nr,1,1])
dY_Vg = np.tile(np.diff(YC, axis=0),[nr,1,1])

xc_dom =XC[plta:pltb, plta:pltb]*1e-3
yc_dom =YC[plta:pltb, plta:pltb]*1e-3

YUg = YC+0.5*dYC # grid where Ug calculated
XUg = XC+0.5*dXC # grid where Vg calculated

Ug_all=np.zeros((nit,nr,ny,nx))
Vg_all=np.zeros((nit,nr,ny,nx))

Ut_all=np.zeros((nit,nr,ny,nx))
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

# To compute W_stern
dyU = np.zeros((ny,nx));
dxV = np.zeros((ny,nx));

zeta = np.zeros((ny,nx));
zeta_intx = np.zeros((ny,nx));
zeta_inty = np.zeros((ny,nx));

termx = np.zeros((ny,nx));
termy = np.zeros((ny,nx));
#
dytermx = np.zeros((ny,nx));
dxtermy = np.zeros((ny,nx));
W_stern = np.zeros((ny,nx));
W_sternall = np.zeros((int(nit/8+1),ny,nx));
Mx_WTall = np.zeros((int(nit/8+1),ny,nx));
My_WTall = np.zeros((int(nit/8+1),ny,nx));
Vort_cndall = np.zeros((int(nit/8+1),ny,nx));
Omega_all = np.zeros((int(nit/8+1),ny,nx));

W5d=np.zeros((nit,nr,ny,nx))

dstart = day_s #5
dend = day_e # (itrs[-1]+1)*ts/720  #10
#itrs= [3600]
idepth = 99 # z indice

# time loop
for it in itrs:
    U = mit.rdmds('U',((it+1)*ts))
    V = mit.rdmds('V',((it+1)*ts))
    W = mit.rdmds('W',((it+1)*ts))
    Vort = vorticity(U[0,:,:],V[0,:,:],dXC,dYC)
    W5d[int((it-itrs[0])), :,:,:]= W;
#
#    P = mit.rdmds('PH',((it+1)*ts))
#
    #dxP = np.diff(P, axis=2)/np.tile(dXC[:,1:],[nr,1,1])
#    dxP = np.diff(P, axis=2)/dX_Ug
#    Vg = dxP/f0
#    dxP300 = np.diff(P[idepth,:,:], axis=1)/np.diff(XC, axis=1)
#    Vg300 = dxP300/f0
#
    #dyP = np.diff(P, axis=1)/np.tile(dYC[1:,:],[nr,1,1])
#    dyP = np.diff(P, axis=1)/dY_Vg
#    Ug = -dyP/f0
#    dyP300 = np.diff(P[idepth,:,:], axis=0)/np.diff(YC, axis=0)
#    Ug300 = -dyP300/f0

# compute U ageostrophic from U bottom
    u_ek = U - U[idepth,:,:]
    v_ek = V - V[idepth,:,:] 
#
#Interpolate Ug in cell center, loosing first grids(bottom and left)
    for z in np.arange(0,nr):
        for y in np.arange(0,ny):
#            Vg_c[z,y,:] = np.interp(XC[y,:],XUg[y,:-1],Vg[z,y,:])
#            Vg300_c[y,:] = np.interp(XC[y,:],XUg[y,:-1],Vg300[y,:])
            Ut_c[z,y,:] = np.interp(XC[y,:],XG[y,:],U[z,y,:])
#            Vort_c[y,:] = np.interp(XC[y,:],XG[y,:],Vort[z,y,:])
        for x in np.arange(0,nx):
#            Ug_c[z,:,x] = np.interp(YC[:,x],YUg[:-1,x],Ug[z,:,x])
#            Ug300_c[:,x] = np.interp(YC[:,x],YUg[:-1,x],Ug300[:,x])
            Vt_c[z,:,x] = np.interp(YC[:,x],YG[:,x],V[z,:,x]) #verify if V points = YUg
#            Vort_c[:,x] = np.interp(YC[:,x],YG[:,x],Vort_c[z,:,x]) #verify if V points = YUg
#    Vort_c = 0.0*Vort
#    for y in np.arange(0,ny):
#        Vort_c[y,:] = np.interp(XC[y,:],XG[y,:],Vort300[y,:])
#    for x in np.arange(0,nx):
#        Vort_c[:,x] = np.interp(YC[:,x],YG[:,x],Vort_c[:,x]) #verify if V points = YUg
#    Vort_cnd = Vort_c/(0.26/(25*1e3)) 

#
#    Vg_all[(it-itrs[0]),:,:,:]= Vg300_c #cell centered Vg
#    Ug_all[(it-itrs[0]),:,:,:]= Ug300_c #cell centered Ug
#    Ut_all[(it-itrs[0]),:,:,:]= Ut_c; #cell centered u_total
#    Vt_all[(it-itrs[0]),:,:,:]= Vt_c; #cell centered v total
#
    u_ekall[(it-itrs[0]),:,:,:]= u_ek; #u-point u_wt
    v_ekall[(it-itrs[0]),:,:,:]= v_ek; #v-point v_wt
#
    if (it % itpd)==0:
        taux = mit.rdmds('diagTAUX',(it*ts))
        tauy = mit.rdmds('diagTAUY',(it*ts))
        """
        for i in range(1,nx):
            for j in range(1,ny):
                dyU[j,i] = (U[1,j,i]-U[1,j-1,i])/dYC[j,i];
                dxV[j,i] = (V[1,j,i]-V[1,j,i-1])/dXC[j,i];
                zeta[j,i] = dxV[j,i]-dyU[j,i];
        """
#
        for i in range(1,nx-1):
            for j in range(1,ny-1):
                zeta_intx[j,i] =(Vort[j,i]+Vort[j+1,i])/2 #vorticity @ u point
                zeta_inty[j,i] =(Vort[j,i]+Vort[j,i+1])/2
#
                termx[j,i] =taux[j,i]/(f0+zeta_intx[j,i])
                termy[j,i] =tauy[j,i]/(f0+zeta_inty[j,i])
#
        for i in range(1,nx-2):
            for j in range(1,ny-2):
                dytermx[(j+1),i] = (termx[j+1,i]-termx[j,i])/dYC[j+1,i];
                dxtermy[j,(i+1)] = (termy[j,i+1]-termy[j,i])/dXC[j,i+1];
                W_stern[(j+1),(i+1)] = (dxtermy[j+1,i+1]-dytermx[j+1,i+1])/rho0;
        W_sternall[int((it-itrs[0])/itpd),:,:]= W_stern; #u-point u_wt
########### 
##### ====  Wenegrat & Thomas ======= ######
        Vort300 = vorticity(U[99,:,:],V[99,:,:],dXC,dYC)
        zetamin = np.min(Vort300)/f0
        zetamax = np.max(Vort300)/f0
        if (abs(zetamin)>zetamax):
            vort_contour = [zetamin*2/3]
        else:
            vort_contour = [zetamax*2/3]
#
# define eddy center
        plt.figure()
        CS = plt.contour(xc_dom,yc_dom,Vort300[plta:pltb,plta:pltb]/f0, vort_contour)
        c =  center_of_mass(CS.allsegs[0][0])
        plt.close()

#transform to polar coord
        thetaC = np.arctan2(YC-c[1]*1e3, XC-c[0]*1e3)
        radC = np.hypot(XC-c[0]*1e3,YC-c[1]*1e3)

# Compute Tangential and Angular velocity
        U_theta = Vt_c*np.cos(thetaC) - Ut_c*np.sin(thetaC)
        Omega = U_theta[99,:,:]/radC/(0.26/(25*1e3))

        Vort_c = 0.0*Vort300
        for y in np.arange(0,ny):
            Vort_c[y,:] = np.interp(XC[y,:],XG[y,:],Vort300[y,:])
        for x in np.arange(0,nx):
            Vort_c[:,x] = np.interp(YC[:,x],YG[:,x],Vort_c[:,x]) #verify if V points = YUg
        Vort_cnd = Vort_c/(0.26/(25*1e3))

# Compute Ekman transport Wenegrat & Thomas
#        Mx_WT = np.max(taux)/(rho0)*Ro*(Vort_cnd-2*Omega)*np.sin(thetaC)*np.cos(thetaC)/((f0+Ro*2*Omega)*(f0+Ro*Vort_cnd)-Ro**2*Omega**2)
#        My_WT = -np.max(taux)/(rho0)*(f0+Ro*Omega+Ro*2*Omega*np.sin(thetaC)**2+Ro*Vort_cnd*np.cos(thetaC)**2)/((f0+Ro*2*Omega)*(f0+Ro*Vort_cnd)-Ro**2*Omega**2)
#        Mx_WT = np.max(taux)*Ro/(rho0*f0)*Ro*(Vort_cnd-2*Omega)*np.sin(thetaC)*np.cos(thetaC)/((1+Ro*2*Omega)*(1+Ro*Vort_cnd)-Ro**2*Omega**2)
#        My_WT = -np.max(taux)/(rho0*f0)*(1+Ro*Omega+Ro*2*Omega*np.sin(thetaC)**2+Ro*Vort_cnd*np.cos(thetaC)**2)/((1+Ro*2*Omega)*(1+Ro*Vort_cnd)-Ro**2*Omega**2)

##        Mx_WT = Ro*(Vort_c-2*Omega)*np.sin(thetaC)*np.cos(thetaC)/((1+Ro*2*Omega)*(1+Ro*Vort_c)-Ro**2*Omega**2)
##        My_WT = -(1+Ro*Omega+Ro*2*Omega*np.sin(thetaC)**2+Ro*Vort_c*np.cos(thetaC)**2)/((1+Ro*2*Omega)*(1+Ro*Vort_c)-Ro**2*Omega**2)
# Small Ro
#        Mx_WT = (np.max(taux)/(rho0*f0))*Ro*(Vort_cnd-2*Omega)*np.sin(thetaC)*np.cos(thetaC)
#        My_WT = -(np.max(taux)/(rho0*f0))*(1-Ro*((Vort_cnd-Omega)*np.sin(thetaC)**2 + Omega*np.cos(thetaC)**2))
#
#        Mx_WTall[int((it-itrs[0])/itpd),:,:]=Mx_WT # @cell center
#        My_WTall[int((it-itrs[0])/itpd),:,:]=My_WT # @cell center
        Vort_cndall[int((it-itrs[0])/itpd),:,:]=Vort_cnd # @cell center
        Omega_all[int((it-itrs[0])/itpd),:,:]=Omega # @cell center


#======= time loop end
# W time averaged
W_ta = np.mean(W5d, axis=0)

Omega_m = np.mean(Omega_all, axis=0)
Vort_m = np.mean(Vort_cndall, axis=0)
Mx_WTm = (np.max(taux)/(rho0*f0))*Ro*(Vort_m-2*Omega_m)*np.sin(thetaC)*np.cos(thetaC)
My_WTm = -(np.max(taux)/(rho0*f0))*(1-Ro*((Vort_m-Omega_m)*np.sin(thetaC)**2 + Omega*np.cos(thetaC)**2))
#
# Compute ekman transport
Mx_ek = np.trapz(u_ekall[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1)
My_ek = np.trapz(v_ekall[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1)

# Compute the divergence
#dyMy_WT = np.diff(My_WTall, axis=1)/np.tile(dYC[1:,:],[6,1,1]) # V-point
dxMx_WT = np.diff(Mx_WTm, axis=1)/dXC[:,1:] # U-point
dyMy_WT = np.diff(My_WTm, axis=0)/dYC[1:,:] # V-point

dxMx_ek = np.diff(Mx_ek, axis=2)/np.tile(dXC[:,1:],[nit,1,1]) # cell center
dyMy_ek = np.diff(My_ek, axis=1)/np.tile(dYC[1:,:],[nit,1,1]) #
"""
#interpolate to cell center
dxMx_WTi=np.zeros((nit,ny,nx))
dyMy_WTi=np.zeros((nit,ny,nx))
for it in np.arange(0,int(nit/itpd)):
    for y in np.arange(0,ny):
#        dxMx_WTi[it,y,:] = np.interp(XC[y,:],XUg[y,:-1],dxMx_WT[it,y,:])
        dxMx_WTi[it,y,:] = np.interp(XC[y,:],XG[y,1:],dxMx_WT[it,y,:])
    for x in np.arange(0,nx):
#        dyMy_WTi[it,:,x] = np.interp(YC[:,x],YUg[:-1,x],dyMy_WT[it,:,x])
        dyMy_WTi[it,:,x] = np.interp(YC[:,x],YG[1:,x],dyMy_WT[it,:,x])
"""
#interpolate to cell center
dxMx_WTi=np.zeros((ny,nx))
dyMy_WTi=np.zeros((ny,nx))
#for it in np.arange(0,int(nit/itpd)):
for y in np.arange(0,ny):
#    dxMx_WTi[it,y,:] = np.interp(XC[y,:],XUg[y,:-1],dxMx_WT[it,y,:])
    dxMx_WTi[y,:] = np.interp(XC[y,:],XG[y,1:],dxMx_WT[y,:])
for x in np.arange(0,nx):
#        dyMy_WTi[it,:,x] = np.interp(YC[:,x],YUg[:-1,x],dyMy_WT[it,:,x])
    dyMy_WTi[:,x] = np.interp(YC[:,x],YG[1:,x],dyMy_WT[:,x])


# Normalization
taux0 = mit.rdmds('diagTAUX',1440)
tau0 = np.max(taux)
Roe = np.max(u_wtall)/(Rmax*f0)

#W_B = (dxMx_i + dyMy_i)
W_ek = dxMx_ek[:,1:,:] + dyMy_ek[:,:,1:] # no need to interpolate because already at cell center
W_WTi = dxMx_WTi + dyMy_WTi


#W_Bm = np.mean(W_B, axis=0)
#Mx_WTm = np.mean(Mx_WTall, axis=0)
#My_WTm = np.mean(My_WTall, axis=0)
#Mx_WTm = Mx_WTall[3,:,:]
#My_WTm = My_WTall[3,:,:]

W_ekman = np.mean(W_ek, axis=0)#/(tau0*Roe/(rho0*f0*Rmax))
Mx_ekman = np.mean(Mx_ek, axis=0)#/(tau0*Roe/(rho0*f0))
My_ekman = np.mean(My_ek, axis=0)#/(tau0/(rho0*f0))

Wm_stern=np.mean(W_sternall, axis=0)
#W_WT = np.mean(W_WTi, axis=0)#/(tau0*Roe/(rho0*f0*Rmax))
#W_WT = W_WTi[3,:,:]#/(tau0*Roe/(rho0*f0*Rmax))

W_WT = W_WTi
###  ============ PLOT ============
# W_ekman bottom 
#
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
#
time = (itrs)*dumpfreq/3600
maxdepthi = np.where(W_ta==np.max(W_ta))
maxdepth = maxdepthi[0][0]
#
#
## ===================  Plot Ekman pumping ===============================
we_range = np.linspace(-5.0, 5.0, 101, endpoint=True)
we_ticks = np.linspace(-5.0, 5.0, 11, endpoint=True)
fig = plt.figure(figsize=(15,8))
fig.suptitle('Ekman pumping, day %d to %d' % (dstart, dend),fontsize=12)
#
ax1 = plt.subplot(231)
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
# ===============================================
ax2 = plt.subplot(232, sharex=ax1)
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
# ===============================================
#
ax3 = plt.subplot(233, sharex=ax1)
ax3.set_aspect(1)
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

# ================================================
#
ax4 = plt.subplot(236, sharex=ax1)
plt.contourf(xc_dom,yc_dom,W_WT[plta:pltb,plta:pltb]*conv_factor,we_range,cmap=cm.seismic)#
#plt.contourf(xc_dom,yc_dom,W_Bm[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(ticks=we_ticks, format='%1.1f')
cb.set_label(r'$W \, (m/day)$', labelpad=-40, y=1.1, rotation=0)
rmax1 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*1e-3, color='r',ls=':', fill=False)
ax3.add_artist(rmax1)
rmax12 = plt.Circle((XC[icy,icx]*1e-3, YC[icy,icx]*1e-3), Rmax*2e-3, color='r',ls=':', fill=False)
ax4.add_artist(rmax12)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
ax4.set_aspect(1)
#plt.text(posext,posext,'$W_{extrema}=$ [%1.1f, %1.1f] $m/day$' % (np.nanmin(W_WT)*conv_factor,np.nanmax(W_WT)*conv_factor), fontsize=10)
#plt.text(posext,posext,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(W_Bm)*1e3,np.max(W_Bm)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
plt.title("Ekman pumping, Wenegrat&Thomas",fontsize=10)
# ===============================================
#
ax5 = plt.subplot(223)
#ax5.set_aspect(1)
plt.plot(yc_dom[:,icx],W_ta[maxdepth,plta:pltb,icx]*conv_factor, label='Time averaged')
plt.plot(yc_dom[:,icx],W_ekman[plta:pltb,icx]*conv_factor, label=r'$\overline{\partial_x \int \, U_e + \partial_y \int \, V_e}$')
plt.plot(yc_dom[:,icx],Wm_stern[plta:pltb,icx]*conv_factor, label='Stern')
plt.plot(yc_dom[:,icx],W_WT[plta:pltb,icx]*conv_factor, label='Wenegrat&Thomas')
plt.legend(loc='best')
plt.grid()
plt.xlabel("y (km)")
plt.ylabel("Vertical velocity (m/day)")
#plt.title("Comparison",fontsize=10)
#
#plt.tight_layout(pad=1)
plt.savefig('./figures/W_ek_mday_day%d-%d.png' % (dstart, dend))

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
plt.title(r'$\overline{M_e^x}=\overline{\int \, U_e}$ day %d-%d' % (dstart, dend),fontsize=11)
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
plt.title(r'$\overline{M_e^y}=\overline{\int \, V_e}$ day %d-%d' % (dstart, dend),fontsize=11)
#plt.text(520,-320,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
plt.tight_layout(pad = 1)
plt.savefig('./figures/Me_day%d-%d.png' % (dstart, dend))
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
plt.title(r'$\overline{M_e^x}=\overline{\int \, U_e}$ day %d-%d' % (dstart, dend),fontsize=11)
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
plt.title(r'$\overline{M_e^y}=\overline{\int \, V_e}$ day %d-%d' % (dstart, dend),fontsize=11)
#plt.text(520,-320,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
plt.tight_layout (pad = 1)
plt.savefig('./figures/Me_WT_day%d-%d.png' % (dstart, dend))

#
#####################################################################
