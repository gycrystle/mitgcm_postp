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
plta = 30
pltb = 300
idepth = 24

timestep = 120 #timestep input in second
dumpfreq = 7200 # in second

day_s = 0
day_e = 1

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

dstart = day_s #5
dend = day_e # (itrs[-1]+1)*ts/720  #10
#itrs= [3600]
#idepth = 99 # z indice for integration

# time loop
for it in itrs:
    U = mit.rdmds('U',((it+1)*ts))
    """
    # Sometimes for some reason 1 time step file cannot be read, 
    # this just a "quick solution" to skip the broken timestep
    if (it+1)*30==1080:
        W = mit.rdmds('W',((it+2)*ts))
    else:
        W = mit.rdmds('W',((it+1)*ts))
    """
#
    V = mit.rdmds('V',((it+1)*ts))
#
    P = mit.rdmds('PH',((it+1)*ts))
#
    #dxP = np.diff(P, axis=2)/np.tile(dXC[:,1:],[nr,1,1])
#    dxP = np.diff(P, axis=2)/dX_Ug
#    Vg = dxP/f0
    dxP300 = np.diff(P[idepth,:,:], axis=1)/np.diff(XC, axis=1)
    Vg300 = dxP300/f0
#
    #dyP = np.diff(P, axis=1)/np.tile(dYC[1:,:],[nr,1,1])
#    dyP = np.diff(P, axis=1)/dY_Vg
#    Ug = -dyP/f0
    dyP300 = np.diff(P[idepth,:,:], axis=0)/np.diff(YC, axis=0)
    Ug300 = -dyP300/f0
# compute U ageostrophic from U bottom
    u_wt = U - U[idepth,:,:]
    v_wt = V - V[idepth,:,:] 
#
#Interpolate Ug in cell center, losing first grids(bottom and left)
    for z in np.arange(0,nr):
        for y in np.arange(0,ny):
#            Vg_c[z,y,:] = np.interp(XC[y,:],XUg[y,:-1],Vg[z,y,:])
            Vg300_c[y,:] = np.interp(XC[y,:],XUg[y,:-1],Vg300[y,:])
            Ut_c[z,y,:] = np.interp(XC[y,:],XG[y,:],U[z,y,:])
        for x in np.arange(0,nx):
#            Ug_c[z,:,x] = np.interp(YC[:,x],YUg[:-1,x],Ug[z,:,x])
            Ug300_c[:,x] = np.interp(YC[:,x],YUg[:-1,x],Ug300[:,x])
            Vt_c[z,:,x] = np.interp(YC[:,x],YG[:,x],V[z,:,x]) #verify if V points = YUg
#
    Vg_all[(it-itrs[0]),:,:,:]= Vg300_c #cell centered Vg
    Ug_all[(it-itrs[0]),:,:,:]= Ug300_c #cell centered Ug
    Ut_all[(it-itrs[0]),:,:,:]= Ut_c; #cell centered u_total
    Vt_all[(it-itrs[0]),:,:,:]= Vt_c; #cell centered v total
#
    u_wtall[(it-itrs[0]),:,:,:]= u_wt; #u-point u_wt
    v_wtall[(it-itrs[0]),:,:,:]= v_wt; #v-point v_wt
#
    if (it % itpd)==0:
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
#
        for i in range(1,nx-2):
            for j in range(1,ny-2):
                dytermx[(j+1),i] = (termx[j+1,i]-termx[j,i])/dYC[j+1,i];
                dxtermy[j,(i+1)] = (termy[j,i+1]-termy[j,i])/dXC[j,i+1];
                W_stern[(j+1),(i+1)] = (dxtermy[j+1,i+1]-dytermx[j+1,i+1])/rho0;
        W_sternall[int((it-itrs[0])/itpd),:,:]= W_stern; #u-point u_wt


#======= time loop end
# Compute cyclogeostrophic velocity
#Ucg = 
# Compute ageostrophic flow
Ue = Ut_all - Ug_all
Ve = Vt_all - Vg_all
# Compute ekman transport
Mx = np.trapz(Ue[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1) 
My = np.trapz(Ve[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1)

Mx_wt = np.trapz(u_wtall[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1)
My_wt = np.trapz(v_wtall[:,0:idepth,:,:],-RC[0:idepth].squeeze(), axis=1)


# Compute the divergence
dxMx = np.diff(Mx, axis=2)/np.tile(dXC[:,1:],[nit,1,1])
dyMy = np.diff(My, axis=1)/np.tile(dYC[1:,:],[nit,1,1])

dxMx_wt = np.diff(Mx_wt, axis=2)/np.tile(dXC[:,1:],[nit,1,1])
dyMy_wt = np.diff(My_wt, axis=1)/np.tile(dYC[1:,:],[nit,1,1])

dxMx_i=np.zeros((nit,ny,nx))
dyMy_i=np.zeros((nit,ny,nx))
for it in np.arange(0,nit):
    for y in np.arange(0,ny):
        dxMx_i[it,y,:] = np.interp(XC[y,:],XUg[y,:-1],dxMx[it,y,:])
    for x in np.arange(0,nx):
        dyMy_i[it,:,x] = np.interp(YC[:,x],YUg[:-1,x],dyMy[it,:,x])

# Normalization
taux0 = mit.rdmds('diagTAUX',1440)
tau0 = np.max(taux)
Roe = np.max(u_wtall)/(Rmax*f0)

W_B = (dxMx_i + dyMy_i)
W_wt = dxMx_wt[:,1:,:] + dyMy_wt[:,:,1:] # no need to interpolate because already at cell center

W_Bm = np.mean(W_B, axis=0)
Mxm = np.mean(Mx, axis=0)
Mym = np.mean(My, axis=0)

Wm_wt = np.mean(W_wt, axis=0)#/(tau0*Roe/(rho0*f0*Rmax))
Mxm_wt = np.mean(Mx_wt, axis=0)#/(tau0*Roe/(rho0*f0))
Mym_wt = np.mean(My_wt, axis=0)#/(tau0/(rho0*f0))

Wm_stern=np.mean(W_sternall, axis=0)

KEg = np.trapz((Ug_all[:,0:idepth,:,:]**2 + Vg_all[:,0:idepth,:,:]**2),axis=1)
KEag = np.trapz((Ue[:,0:idepth,:,:]**2 + Ve[:,0:idepth,:,:]**2),axis=1)
KEwt = np.trapz((u_wtall[:,0:idepth,:,:]**2 + v_wtall[:,0:idepth,:,:]**2),axis=1)

#============ PLOT ==========================
# W_ekman bottom 

wekmax=np.max(W_Bm)*1e3
wekmin=np.min(W_Bm)*1e3
wavminmax=max(abs(wekmin), abs(wekmax))
w_range = np.linspace(-wavminmax, wavminmax, 101, endpoint=True)
wtmax=np.max(Wm_wt)*1e3
wtmin=np.min(Wm_wt)*1e3
wtminmax=max(abs(wtmin), abs(wtmax))
#w_range = np.linspace(-wavminmax, wavminmax, 101, endpoint=True)
wt_range = np.linspace(-wtminmax, wtminmax, 101, endpoint=True)
mminmax = max(abs(np.min(Mxm)), abs(np.max(Mxm)),abs(np.min(Mym)), abs(np.max(Mym)))

mwtminmax = max(abs(np.min(Mxm_wt)), abs(np.max(Mxm_wt)),abs(np.min(Mym_wt)), abs(np.max(Mym_wt)))
m_range = np.linspace(-mminmax, mminmax, 101, endpoint=True)
mwt_range = np.linspace(-mwtminmax, mwtminmax, 101, endpoint=True)
#norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
#
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
#
#time = (itrs+1)*120*ts/3600
time = (itrs)*dumpfreq/3600

#
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

#
# Ekman tranport plot
we_range = np.linspace(-0.1, 0.1, 101, endpoint=True)
fig = plt.figure(figsize=(5,8))
ax1 = plt.subplot(211)
#plt.contourf(xc_dom,yc_dom,Wm_wt[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
plt.contourf(xc_dom,yc_dom,W_Bm[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
ax1.set_aspect(1)
#plt.text(370,370,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(Wm_wt)*1e3,np.max(Wm_wt)*1e3), fontsize=10)
plt.text(370,370,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(W_Bm)*1e3,np.max(W_Bm)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
plt.title(r'$\overline{W_{e}}=\overline{\partial_x \int \, U_e + \partial_y \int \, V_e}$ day %d-%d' % (dstart, dend),fontsize=11)
#
ax2 = plt.subplot(212, sharex=ax1)
    #ax2 = fig.add_subplot(2, 1, 2)
ax2.set_aspect(1)
plt.contourf(xc_dom,yc_dom,Wm_stern[plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
#CS1 = plt.contour(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3, levels, colors='0.6')
#plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)
plt.text(370,370,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(Wm_stern)*1e3,np.max(Wm_stern)*1e3), fontsize=10)
#    plt.text(375,-1250,'$W_{max}=$ %1.3f $mm/s$' % (wmax))
#    plt.text(375,-1350,'$W_{min}=$ %1.3f $mm/s$' % (wmin))

plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
plt.title(r'$\overline{W_{stern}}=\overline{\frac{1}{\rho_0} \nabla \times \left[\frac{\tau}{f+\zeta}\right]}$', fontsize=14)
#plt.text(520,-320,'$W_{extrema}=$[%1.3f, %1.3f] $mm/s$' % (wmin, wmax), fontsize=10)
plt.tight_layout (pad = 1)
plt.savefig('./figures/W_ek_day%d-%d.png' % (dstart, dend))
#    plt.close()



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

plot_perturb = 1
if plot_perturb == 1:
    u_tilde = 0.0*np.ones((nr,ny,nx,nit));
    v_tilde = 0.0*np.ones((nr,ny,nx,nit));
    dyU = 0.0*np.ones((ny,nx,nit));
    dxV = 0.0*np.ones((ny,nx,nit));
    zeta = 0.0*np.ones((ny,nx,nit));

    for it in itrs:
        u_tilde[:,:,:,it] = Uall[:,:,:,it]-Uall[:,:,:,0]
        v_tilde[:,:,:,it] = Vall[:,:,:,it]-Vall[:,:,:,0]
#
        for i in range(1,250):
            for j in range(1,250):
                dyU[j,i,it] = (u_tilde[1,j,i,it]-u_tilde[1,j-1,i,it])/dYC[j,i];
                dxV[j,i,it] = (v_tilde[1,j,i,it]-v_tilde[1,j,i-1,it])/dXC[j,i];
                zeta[j,i,it] = dxV[j,i,it]-dyU[j,i,it];

# Compute Kinetic Energy of perturbation
K=u_tilde**2+v_tilde**2

lnsqrtK=np.log(np.sqrt(K))

#######################Plot the hovmoller diagrams##########################################

## Hovmoller time vs y 
#
time = itrs+1
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
#w_range = np.linspace(-0.6,0.6,61, endpoint=True)
#w_range2 = np.linspace(-0.3,0.3,61, endpoint=True)
w_range = 101
hstart1 = 0 
hend1 = 60

plt.figure(figsize=(12,6))
plt.plot(time,lnsqrtK[0,125,125,:])
plt.ylabel(r'$\ln K^{1/2}$')
plt.xlabel("time (hour)")
plt.title(r'$\ln K^{1/2}$ at $x = 150km, y = 150km$, day %d-%d' % (dstart, dend))
plt.savefig('./figures/K.png')

fig = plt.figure(figsize=(15,6))
#
ax1 = fig.add_subplot(1, 2, 1)
plt.contourf(time[hstart1:hend1].squeeze()*2,YC[:,125]*1e-3,Wall[33,:,125,hstart1:hend1]*1e3,w_range,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
plt.colorbar(label='$W \ [mm/s]$', format='%1.3f')
#CS1 = plt.contour(time[hstart1:hend1].squeeze()*2,YC[:,125]*1e-3,Wall[33,:,125,hstart1:hend1]*1e3, levels, colors='0.6')
#plt.clabel(CS1, fmt='%2.3f', colors='k', fontsize=10)
plt.ylabel("y (km)")
plt.xlabel("time (hour)")
plt.title(r'$W$ at $x = 150km, z\sim 207m$, day %d-%d' % (int(hstart1/12), int(hend1/12)))
"""
#
"""
ax2 = fig.add_subplot(1, 2, 2)
hstart2 = 120
hend2 = 240
plt.contourf(time[hstart2:hend2].squeeze(),YC[:,125]*1e-3,Wall[33,:,125,hstart2:hend2]*1e3,101,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
plt.colorbar(label='$W \ [mm/s]$', format='%1.3f')
CS3 = plt.contour(time[hstart2:hend2].squeeze(),YC[:,125]*1e-3,Wall[33,:,125,hstart2:hend2]*1e3, levels, colors='0.6')
plt.clabel(CS3, fmt='%2.3f', colors='k', fontsize=10)
plt.ylabel("y (km)")
plt.xlabel("time (hour)")
plt.title(r'$W$ at $x = 150km, z\sim 207m$, day %d-%d' % (int(hstart2/24), int(hend2/24)))
"""
#
"""
plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_W207_day%d_%d.png' % (dstart, dend))

#### Hovmoller surface  vorticity
fig = plt.figure(figsize=(12,6))
#
plt.contourf(time.squeeze(),YC[:,125]*1e-3,zeta[:,125,:]*1e3,101,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
plt.colorbar(label='$W \ [mm/s]$', format='%1.3f')

CS2 = plt.contour(time.squeeze(),YC[:,125]*1e-3,zeta[:,125,:]*1e3, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)

plt.ylabel("y (km)")
plt.xlabel("time (hour)")
plt.title('Surface vorticity ($\zeta$) at $x=150km$, day %d-%d' % (dstart, dend))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_Vorticity_day%d_%d.png' % (dstart, dend))

#### Hovmoller surface u tilde
fig = plt.figure(figsize=(12,6))
#
plt.contourf(time.squeeze(),YC[:,125]*1e-3,u_tilde[0,:,125,:]*1e3,101,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
plt.colorbar(label='$W \ [mm/s]$', format='%1.3f')

CS2 = plt.contour(time.squeeze(),YC[:,125]*1e-3,u_tilde[0,:,125,:]*1e3, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)

plt.ylabel("y (km)")
plt.xlabel("time (hour)")
plt.title('U perturbation ($\zeta$) at $x=150km$, day %d-%d' % (dstart, dend))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_utilde_day%d_%d.png' % (dstart, dend))

#### Hovmoller time vs depth
fig = plt.figure(figsize=(12,6))
#
plt.contourf(time.squeeze(),RC.squeeze(),Wall[:,125,125,:]*1e3,101,norm=MidpointNormalize(midpoint=0.),cmap=cm.seismic)
plt.colorbar(label='$W \ [mm/s]$', format='%1.3f')

CS2 = plt.contour(time.squeeze(),RC.squeeze(),Wall[:,125,125,:]*1e3, levels, colors='0.6')
plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)

plt.ylabel("depth (m)")
plt.xlabel("time (hour)")
plt.title('Vertical velocity ($W$) at $x=150km,\  y=150km$, day %d-%d' % (dstart, dend))

plt.tight_layout(pad=1)
plt.savefig('./figures/hovmoller_Wcore_day%d_%d.png' % (dstart, dend))
"""
#
#####################################################################
