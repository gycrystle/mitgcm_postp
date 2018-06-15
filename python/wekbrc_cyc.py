"""
Plot the w_ekman from ek transport divergence
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors

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

# select plot domain
plta = 40
pltb = 290
idepth = 24

# choose how to compute ageostrophic component
ubottom = 0 #ubottom = 1, similar to wenegrat  thomas

timestep = 120 #timestep input in second
dumpfreq = 7200 # in second

day_s = 15
day_e = 20

ts = dumpfreq/timestep  # time step = ts*dt (in second); = 7200 = dumpfreq

nstart = int(day_s*86400/dumpfreq) # integer of dumpfrec, either start from zero or 
nend = int(day_e*86400/dumpfreq) # no of time step 
itpd = int(86400/dumpfreq)

f0 = 1e-4
rho0 = 999.8
Rmax = 25e3
vmax = 0.4

#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(nstart,nend,1)
nit=itrs.size

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')
icx = int(XC[0,:].size/2)
icy = int(YC[:,0].size/2)
#transform to polar coord
thetaC = np.arctan2(YC-YC[icy, icx], XC-XC[icy, icx])
radC = np.hypot(XC-XC[icy, icx],YC-YC[icy, icx])

nr = RC.size
nx = XC[0,:].size
ny = YC[:,0].size

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
Vt_all=np.zeros((nit,nr,ny,nx))

u_wtall=np.zeros((nit,nr,ny,nx))
v_wtall=np.zeros((nit,nr,ny,nx))

Ug_c=np.zeros((nr,ny,nx))
Vg_c=np.zeros((nr,ny,nx))

Ut_c=np.zeros((nr,ny,nx))
Vt_c=np.zeros((nr,ny,nx))
#Ur_cygeo = np.zeros((nr, ny, nx))
U_cg_all= np.zeros((nit,nr,ny,nx))
V_cg_all= np.zeros((nit,nr,ny,nx))

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
W_sternall = np.zeros((int(nit/itpd+1),ny,nx));

termx2 = np.zeros((ny,nx));
termy2 = np.zeros((ny,nx));
#
dytermx2 = np.zeros((ny,nx));
dxtermy2 = np.zeros((ny,nx));
W_stern2 = np.zeros((ny,nx));
W_sternall2 = np.zeros((int(nit/itpd+1),ny,nx));

dstart = day_s #5
dend = day_e # (itrs[-1]+1)*ts/720  #10
#itrs= [3600]
#idepth = 99 # z indice for integration
U0 = mit.rdmds('U',0)
V0 = mit.rdmds('V',0)

# time loop
for it in itrs:
    U = mit.rdmds('U',((it+1)*ts))
#
    V = mit.rdmds('V',((it+1)*ts))
#
    P = mit.rdmds('PH',((it+1)*ts))
#
    #dxP = np.diff(P, axis=2)/np.tile(dXC[:,1:],[nr,1,1])
    dxP = np.diff(P, axis=2)/dX_Ug
    Vg = dxP/f0
#
    #dyP = np.diff(P, axis=1)/np.tile(dYC[1:,:],[nr,1,1])
    dyP = np.diff(P, axis=1)/dY_Vg
    Ug = -dyP/f0
#
# compute U ageostrophic from U bottom or U_0
    if ubottom == 1:
        u_wt = U - U[idepth,:,:]
        v_wt = V - V[idepth,:,:] 
    else:
        u_wt = U - U0
        v_wt = V - V0

#
#Interpolate Ug in cell center, losing first grids(bottom and left)
    for z in np.arange(0,nr):
        for y in np.arange(0,ny):
            Vg_c[z,y,:] = np.interp(XC[y,:],XUg[y,:-1],Vg[z,y,:])
            Ut_c[z,y,:] = np.interp(XC[y,:],XG[y,:],U[z,y,:])
        for x in np.arange(0,nx):
            Ug_c[z,:,x] = np.interp(YC[:,x],YUg[:-1,x],Ug[z,:,x])
            Vt_c[z,:,x] = np.interp(YC[:,x],YG[:,x],V[z,:,x]) #verify if V points = YUg
#
# transform u_geostrophic to cylindrical to compute cyclogeo vel
    U_theta_g = Vg_c*np.cos(thetaC)-Ug_c*np.sin(thetaC)
#    U_r_g = Ug_c*np.cos(thetaC) + Vg_c*np.sin(thetaC)
# Tangential vel component of gradient wind, cell centered
    U_theta_cygeo = (f0*radC/2)-np.sqrt(f0**2*radC**2/4 - f0*radC*U_theta_g)
# transform U_cyclogeo to cartesian
    U_cygeo = - U_theta_cygeo*np.sin(thetaC) # Ur_cygeo*np.cos(thetaC)
    V_cygeo = U_theta_cygeo*np.cos(thetaC) # Ur_cygeo*np.sin(thetaC)
#
    Vg_all[(it-itrs[0]),:,:,:]= Vg_c#np.tile(Vg300_c,(nr,1,1)) #cell centered Vg
    Ug_all[(it-itrs[0]),:,:,:]= Ug_c#np.tile(Ug300_c,(nr,1,1)) #cell centered Ug
    Ut_all[(it-itrs[0]),:,:,:]= Ut_c; #cell centered u_total
    Vt_all[(it-itrs[0]),:,:,:]= Vt_c; #cell centered v total
    U_cg_all[(it-itrs[0]),:,:,:]= U_cygeo #np.tile(U_cygeo, (nr,1,1));
    V_cg_all[(it-itrs[0]),:,:,:]= V_cygeo #np.tile(V_cygeo, (nr,1,1));
#
    u_wtall[(it-itrs[0]),:,:,:]= u_wt; #u-point u_wt
    v_wtall[(it-itrs[0]),:,:,:]= v_wt; #v-point v_wt
#
    if it != 0 and (it % itpd)==0:
        taux = mit.rdmds('diagTAUX',(it*ts))
        tauy = mit.rdmds('diagTAUY',(it*ts))

        for i in range(1,nx):
            for j in range(1,ny):
                dyU[j,i] = (U[1,j,i]-U[1,j-1,i])/dYC[j,i];
                dxV[j,i] = (V[1,j,i]-V[1,j,i-1])/dXC[j,i];
                zeta[j,i] = dxV[j,i]-dyU[j,i];
#
        for i in range(1,nx-1):
            for j in range(1,ny-1):
                zeta_intx[j,i] =(zeta[j,i]+zeta[j+1,i])/2
                zeta_inty[j,i] =(zeta[j,i]+zeta[j,i+1])/2
#
                termx[j,i] =taux[j,i]/(1e-4+zeta_intx[j,i])
                termy[j,i] =tauy[j,i]/(1e-4+zeta_inty[j,i])
                termx2[j,i] =taux[j,i]/(1e-4+2*zeta_intx[j,i])
                termy2[j,i] =tauy[j,i]/(1e-4+2*zeta_inty[j,i])
#
        for i in range(1,nx-2):
            for j in range(1,ny-2):
                dytermx[(j+1),i] = (termx[j+1,i]-termx[j,i])/dYC[j+1,i];
                dxtermy[j,(i+1)] = (termy[j,i+1]-termy[j,i])/dXC[j,i+1];
                dytermx2[(j+1),i] = (termx2[j+1,i]-termx2[j,i])/dYC[j+1,i];
                dxtermy2[j,(i+1)] = (termy2[j,i+1]-termy2[j,i])/dXC[j,i+1];
#
                W_stern[(j+1),(i+1)] = (dxtermy[j+1,i+1]-dytermx[j+1,i+1])/rho0;
                W_stern2[(j+1),(i+1)] = (dxtermy2[j+1,i+1]-dytermx2[j+1,i+1])/rho0;
                
        W_sternall[int((it-itrs[0])/itpd),:,:]= W_stern; #u-point u_wt
        W_sternall2[int((it-itrs[0])/itpd),:,:]= W_stern2; #


#======= time loop end

################## Compute ageostrophic flow ##############

#1 from geostrophic balance
Ue_g = Ut_all - Ug_all
Ve_g = Ut_all - Vg_all


#2 from cyclogeostrophic balance
Ue_cyg = Ut_all - U_cg_all
Ve_cyg = Vt_all - V_cg_all

# Compute ekman transport

Mx_g = np.trapz(Ue_g[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1) 
My_g = np.trapz(Ve_g[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1)

Mx_cyg = np.trapz(Ue_cyg[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1)
My_cyg = np.trapz(Ve_cyg[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1)

Mx_wt = np.trapz(u_wtall[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1)
My_wt = np.trapz(v_wtall[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1)

# Compute the divergence
dxMx_g = np.diff(Mx_g, axis=2)/np.tile(dXC[:,1:],[nit,1,1])
dyMy_g = np.diff(My_g, axis=1)/np.tile(dYC[1:,:],[nit,1,1])

dxMx_cyg = np.diff(Mx_cyg, axis=2)/np.tile(dXC[:,1:],[nit,1,1])
dyMy_cyg = np.diff(My_cyg, axis=1)/np.tile(dYC[1:,:],[nit,1,1])

dxMx_wt = np.diff(Mx_wt, axis=2)/np.tile(dXC[:,1:],[nit,1,1])
dyMy_wt = np.diff(My_wt, axis=1)/np.tile(dYC[1:,:],[nit,1,1])

dxMx_gi=np.zeros((nit,ny,nx))
dyMy_gi=np.zeros((nit,ny,nx))

dxMx_cygi=np.zeros((nit,ny,nx))
dyMy_cygi=np.zeros((nit,ny,nx))

for it in np.arange(0,nit):
    for y in np.arange(0,ny):
        dxMx_gi[it,y,:] = np.interp(XC[y,:],XUg[y,:-1],dxMx_g[it,y,:])
        dxMx_cygi[it,y,:] = np.interp(XC[y,:],XUg[y,:-1],dxMx_cyg[it,y,:])
    for x in np.arange(0,nx):
        dyMy_gi[it,:,x] = np.interp(YC[:,x],YUg[:-1,x],dyMy_g[it,:,x])
        dyMy_cygi[it,:,x] = np.interp(YC[:,x],YUg[:-1,x],dyMy_cyg[it,:,x])

# Normalization coeff
"""
taux0 = mit.rdmds('diagTAUX',1440)
tau0 = np.max(taux)
Roe = np.max(u_wtall)/(Rmax*f0)
"""

# compute time averaged Ekman pumping
We_g = np.mean((dxMx_gi + dyMy_gi), axis=0)
We_cyg = np.mean((dxMx_cygi + dyMy_cygi), axis=0)
We_wt = np.mean((dxMx_wt[:,1:,:] + dyMy_wt[:,:,1:]), axis=0) # no need to interpolate because already at cell center
We_stern=np.mean(W_sternall, axis=0)
We_stern2=np.mean(W_sternall2, axis=0)

#W_Bm = np.mean(W_B, axis=0)

# time averaging ekman transports

Mx_g = np.mean(Mx_g, axis=0)
My_g = np.mean(My_g, axis=0)

Mx_g300 = np.mean(Mx_g300, axis=0)
My_g300 = np.mean(My_g300, axis=0)

Mx_cyg = np.mean(Mx_cyg, axis=0)
My_cyg = np.mean(My_cyg, axis=0)

#Wm_wt = np.mean(W_wt, axis=0)#/(tau0*Roe/(rho0*f0*Rmax))
Mx_wt = np.mean(Mx_wt, axis=0)#/(tau0*Roe/(rho0*f0))
My_wt = np.mean(My_wt, axis=0)#/(tau0/(rho0*f0))

"""
KEg = np.trapz((Ug_all[:,0:idepth,:,:]**2 + Vg_all[:,0:idepth,:,:]**2),axis=1)
KEag = np.trapz((Ue[:,0:idepth,:,:]**2 + Ve[:,0:idepth,:,:]**2),axis=1)
KEwt = np.trapz((u_wtall[:,0:idepth,:,:]**2 + v_wtall[:,0:idepth,:,:]**2),axis=1)
"""
#============ PLOT ==========================
# Ekman pumping 

wekmax=np.max(We_g)*1e3
wekmin=np.min(We_g)*1e3
wavminmax=max(abs(wekmin), abs(wekmax))
w_range = np.linspace(-wavminmax, wavminmax, 101, endpoint=True)
wtmax=np.max(We_wt)*1e3
wtmin=np.min(We_wt)*1e3
wtminmax=max(abs(wtmin), abs(wtmax))
#w_range = np.linspace(-wavminmax, wavminmax, 101, endpoint=True)
#wt_range = np.linspace(-wtminmax, wtminmax, 101, endpoint=True)
#mminmax = max(abs(np.min(Mxm)), abs(np.max(Mxm)),abs(np.min(Mym)), abs(np.max(Mym)))

#mwtminmax = max(abs(np.min(Mxm_wt)), abs(np.max(Mxm_wt)),abs(np.min(Mym_wt)), abs(np.max(Mym_wt)))
#m_range = np.linspace(-mminmax, mminmax, 101, endpoint=True)
#mwt_range = np.linspace(-mwtminmax, mwtminmax, 101, endpoint=True)
#norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
#
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
#
#time = (itrs+1)*120*ts/3600
time = (itrs)*dumpfreq/3600

#
"""
fig = plt.figure(figsize=(15,6.5))
ax1 = fig.add_subplot(1, 2, 1)
plt.contourf(xc_dom,yc_dom,W_Bm[plta:pltb,plta:pltb]*1e3,w_range,cmap=cm.seismic)
plt.colorbar(label='$\overline{W} \ [mm/s]$', format='%1.3f')
plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (wekmin, wekmax), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$U_e = U_t - U_g$ where $U_g$ calculated from geostrophic balance', fontsize=11)
#
ax2 = fig.add_subplot(1, 2, 2)
plt.contourf(xc_dom,yc_dom,Wm_wt[plta:pltb,plta:pltb]*1e3,wt_range,cmap=cm.seismic)
plt.colorbar(label='$\overline{W} \ [mm/s]$', format='%1.3f')
plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (wtmin, wtmax), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$U_e$ calculated from u - u(z=300)', fontsize=11)
#
plt.suptitle('5 days averaged Ekman pumping $\overline{W_e}$, day %d-%d' % (dstart, dend), fontsize=14)
plt.tight_layout(pad=1)
plt.subplots_adjust(top=0.89)
plt.savefig('./figures/W_ekman%d_%d.png' % (dstart, dend))
#============================================================
#Ekman transport Mx
fig  = plt.figure(figsize=(15,6.5))
ax1 = fig.add_subplot(1, 2, 1)
plt.contourf(xc_dom,yc_dom,Mxm[plta:pltb, plta:pltb],m_range,cmap=cm.seismic)
plt.colorbar(label='$\overline{M_x} \ [m^2/s]$', format='%1.3f')
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$U_e = U_t - U_g$ where $U_g$ calculated from geostrophic balance', fontsize=11)
#
ax2 = fig.add_subplot(1, 2, 2)
plt.contourf(xc_dom,yc_dom,Mxm_wt[plta:pltb, plta:pltb],100, norm=MidpointNormalize(midpoint=0.0),cmap=cm.seismic)
plt.colorbar(label='$\overline{M_x} \ [m^2/s]$', format='%1.3f')
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$U_e$ calculated from u - u(z=300)', fontsize=11)
#
plt.suptitle('5 days averaged Ekman transport $\overline{M_x}$, day %d-%d' % (dstart, dend), fontsize=14)
plt.tight_layout(pad=1)
plt.subplots_adjust(top=0.89)
plt.savefig('./figures/Mx_day%d_%d.png' % (dstart, dend))

#Ekman transport My
fig  = plt.figure(figsize=(15,6.5))
ax1 = fig.add_subplot(1, 2, 1)
plt.contourf(xc_dom,yc_dom,Mym[plta:pltb, plta:pltb],m_range,cmap=cm.seismic)
plt.colorbar(label='$\overline{M_x} \ [m^2/s]$', format='%1.3f')
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$V_e = V_t - V_g$ where $V_g$ calculated from geostrophic balance', fontsize=11)
#
ax2 = fig.add_subplot(1, 2, 2)
plt.contourf(xc_dom,yc_dom,Mym_wt[plta:pltb, plta:pltb],101,norm=MidpointNormalize(midpoint=Mym_wt[plta,plta]),cmap=cm.seismic)
plt.colorbar(label='$\overline{M_x} \ [m^2/s]$', format='%1.3f')
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$V_e$ calculated from v - v(z=300)', fontsize=11)
#
plt.suptitle('5 days averaged Ekman transport $\overline{M_y}$, day %d-%d' % (dstart, dend), fontsize=14)
plt.tight_layout(pad=1)
plt.subplots_adjust(top=0.89)
plt.savefig('./figures/My_day%d_%d.png' % (dstart, dend))

# Kinetic Energy 
plt.figure(figsize=(12,6))
plt.plot(time,KEg[:,int(nx/2+20),int(nx/2+20)], label='geostrophic')
#plt.plot(time,KEag[:,int(nx/2+20),int(nx/2+20)], label='ageostrophic')
plt.plot(time,KEwt[:,int(nx/2+20),int(nx/2+20)], label='ageostrophic -u(300)')
plt.ylabel(r'$\propto KE$')
plt.xlabel("time (hour)")
plt.legend(loc='best', fontsize=10)
plt.title(r'$\propto KE$ at 25km from eddy center, day %d-%d' % (dstart, dend))
plt.savefig('./figures/KE.png')
"""
######################################################################################
# Ekman tranport plot
we_range = np.linspace(-0.1, 0.1, 101, endpoint=True)
fig = plt.figure(figsize=(15,8))
ax1 = plt.subplot(231)
#plt.contourf(xc_dom,yc_dom,Wm_wt[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
plt.contourf(xc_dom,yc_dom,We_g[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
ax1.set_aspect(1)
plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(We_g)*1e3,np.max(We_g)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$\overline{W_{e}}=\overline{\partial_x \int \, U_eg + \partial_y \int \, V_eg}$ day %d-%d' % (dstart, dend),fontsize=11)
#
ax2 = plt.subplot(232, sharex=ax1)
ax2.set_aspect(1)
plt.contourf(xc_dom,yc_dom,We_g300[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(We_g300)*1e3,np.max(We_g300)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
plt.title(r'$\overline{W_{e}}=\overline{\partial_x \int \, U_eg300 + \partial_y \int \, V_eg300}$ day %d-%d' % (dstart, dend),fontsize=11)

ax3 = plt.subplot(233, sharex=ax1)
ax3.set_aspect(1)
plt.contourf(xc_dom,yc_dom,We_wt[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(We_wt)*1e3,np.max(We_wt)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$\overline{W_{e}}=\overline{\partial_x \int \, U_e300 + \partial_y \int \, V_e300}$ day %d-%d' % (dstart, dend),fontsize=11)

ax4 = plt.subplot(234, sharex=ax1)
ax4.set_aspect(1)
plt.contourf(xc_dom,yc_dom,We_cyg[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(We_cyg)*1e3,np.max(We_cyg)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$\overline{W_{e}}=\overline{\partial_x \int \, U_{ecg} + \partial_y \int \, V_{ecg}}$ day %d-%d' % (dstart, dend),fontsize=11)

ax5 = plt.subplot(235, sharex=ax1)
ax5.set_aspect(1)
plt.contourf(xc_dom,yc_dom,We_stern[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
#CS1 = plt.contour(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3, levels, colors='0.6')
#plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)
plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(We_stern)*1e3,np.max(We_stern)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
plt.title(r'$\overline{W_{stern}}=\overline{\frac{1}{\rho_0} \nabla \times \left[\frac{\tau}{f+\zeta}\right]}$', fontsize=14)
###
ax6 = plt.subplot(236, sharex=ax1)
ax6.set_aspect(1)
plt.contourf(xc_dom,yc_dom,We_stern2[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(We_stern2)*1e3,np.max(We_stern2)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$\overline{W_{stern*}}=\overline{\frac{1}{\rho_0} \nabla \times \left[\frac{\tau}{f+2\zeta}\right]}$', fontsize=14)

plt.tight_layout (pad = 1)
plt.savefig('./figures/W_ek_day%d-%d.png' % (dstart, dend))
#    plt.close()


"""
#Ekman Transport plot
#we_range = np.linspace(-0.02, 0.02, 101, endpoint=True)
fig = plt.figure(figsize=(5,8))
ax1 = plt.subplot(211)
plt.contourf(xc_dom,yc_dom,Mxm[plta:pltb, plta:pltb],101,norm=MidpointNormalize(midpoint=Mxm[plta,plta]),cmap=cm.seismic)
cb = plt.colorbar(format='%1.3f')
cb.set_label(r'$M_e \, (m^2/s)$', labelpad=-40, y=1.1, rotation=0)
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
    #ax2 = fig.add_subplot(2, 1, 2)
ax2.set_aspect(1)
plt.contourf(xc_dom,yc_dom,Mym[plta:pltb, plta:pltb],101,norm=MidpointNormalize(midpoint=Mym[plta,plta]),cmap=cm.seismic)
cb=plt.colorbar(format='%1.3f')
cb.set_label(r'$M_e \, (m^2/s)$', labelpad=-40, y=1.1, rotation=0)
#CS1 = plt.contour(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3, levels, colors='0.6')
#plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)
#plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(Wm_stern)*1e3,np.max(Wm_stern)*1e3), fontsize=10)
#    plt.text(375,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
#    plt.text(375,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
plt.title(r'$\overline{M_e^y}=\overline{\int \, V_e}$ day %d-%d' % (dstart, dend),fontsize=11)
#plt.text(520,-320,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
plt.tight_layout (pad = 1)
plt.savefig('./figures/Me_day%d-%d.png' % (dstart, dend))
#    plt.close()
"""

#####################################################################
