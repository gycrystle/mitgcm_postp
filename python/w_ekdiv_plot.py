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

ts = 60  # time step = ts*dt (in second); = 7200 = dumpfreq
nr = 100  # no of grid vertical 
nx = 180 # no of grid in x
ny = 180 # no of grid in y
endtime = 3600 # in time step = 120s
nstart = 0 # in day*2/24
nend = int(endtime/ts) # no of time step 
f0 = 1e-4

# select plot domain
plta = 31
pltb = 151
xc_dom =XC[plta:pltb, plta:pltb]*1e-3
yc_dom =YC[plta:pltb, plta:pltb]*1e-3

#define time step to plot (one time step is 120s, 1 day = 86400 s = 720 time step
itrs=np.arange(nstart,nend,1)
nit=itrs.size

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')
RC = mit.rdmds('RC')
XG = mit.rdmds('XG')
YG = mit.rdmds('YG')
dXC = mit.rdmds('DXC')
dYC = mit.rdmds('DYC')
dXG = mit.rdmds('DXG')
dYG = mit.rdmds('DYG')
dX_Ug = np.tile(np.diff(XC, axis=1),[nr,1,1])
dY_Vg = np.tile(np.diff(YC, axis=0),[nr,1,1])

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

dstart = (itrs[0]+1)*ts/720 #5
dend = (itrs[-1]+1)*ts/720  #10
#itrs= [3600]
idepth = 99 # z indice for integration

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
    dxP = np.diff(P, axis=2)/dX_Ug
    Vg = dxP/f0
#
    #dyP = np.diff(P, axis=1)/np.tile(dYC[1:,:],[nr,1,1])
    dyP = np.diff(P, axis=1)/dY_Vg
    Ug = -dyP/f0

# compute U ageostrophic from U bottom
    u_wt = U - U[99,:,:]
    v_wt = V - V[99,:,:] 
#
#Interpolate Ug in cell center, losing first grids(bottom and left)
    for z in np.arange(0,nr):
        for y in np.arange(0,ny):
            Vg_c[z,y,:] = np.interp(XC[y,:],XUg[y,:-1],Vg[z,y,:])
            Ut_c[z,y,:] = np.interp(XC[y,:],XG[y,:],U[z,y,:])
        for x in np.arange(0,nx):
            Ug_c[z,:,x] = np.interp(YC[:,x],YUg[:-1,x],Ug[z,:,x])
            Vt_c[z,:,x] = np.interp(YC[:,x],YG[:,x],V[z,:,x]) #verify if V points = YUg
    Vg_all[(it-itrs[0]),:,:,:]= Vg_c #cell centered Vg
    Ug_all[(it-itrs[0]),:,:,:]= Ug_c #cell centered Ug
    Ut_all[(it-itrs[0]),:,:,:]= Ut_c; #cell centered u_total
    Vt_all[(it-itrs[0]),:,:,:]= Vt_c; #cell centered v total
#
    u_wtall[(it-itrs[0]),:,:,:]= u_wt; #u-point u_wt
    v_wtall[(it-itrs[0]),:,:,:]= v_wt; #v-point v_wt

#======= time loop end

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
        dxMx_i[it,:,x] = np.interp(YC[:,x],YUg[:-1,x],dyMy[it,:,x])

#W_B = dxMx[:,1:,:] + dyMy[:,:,1:]
W_B = dxMx_i + dyMy_i
W_wt = dxMx_wt + dyMy_wt # no need to interpolate because already at cell center

W_Bm = np.mean(W_B, axis=0)
Mxm = np.mean(Mx, axis=0)
Mym = np.mean(My, axis=0)

Wm_wt = np.mean(W_wt, axis=0)
Mxm_wt = np.mean(Mx_wt, axis=0)
Mym_wt = np.mean(My_wt, axis=0)

KEg = np.trapz((Ug_all[:,0:idepth,:,:]**2 + Vg_all[:,0:idepth,:,:]**2),axis=1)
KEag = np.trapz((Ue[:,0:idepth,:,:]**2 + Ve[:,0:idepth,:,:]**2),axis=1)

#============ PLOT ==========================
# W_ekman bottom 

wekmax=np.max(W_Bm)#*1e3
wekmin=np.min(W_Bm)#*1e3
wavminmax=max(abs(wekmin), abs(wekmax))
w_range = np.linspace(-wavminmax, wavminmax, 101, endpoint=True)
#w_range = np.linspace(-wavminmax, wavminmax, 101, endpoint=True)
wt_range = 100
#
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
#
time = (itrs+1)*120*ts/3600
#
fig = plt.figure(figsize=(14,6))
ax1 = fig.add_subplot(1, 2, 1)
plt.contourf(xc_dom,yc_dom,W_Bm[plta:pltb,plta:pltb]*1e3,w_range,cmap=cm.seismic)
plt.colorbar(label='$\overline{W} \ [mm/s]$', format='%1.3f')
plt.text(50,120,'$W_{max}=$ %1.3f $mm/s$' % (wekmax))
plt.text(50,50,'$W_{min}=$ %1.3f $mm/s$' % (wekmin))
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$U_e = U_t - U_g$ where $U_g$ calculated from geostrophic balance')
#
ax2 = fig.add_subplot(1, 2, 2)
plt.contourf(xc_dom,yc_dom,Wm_wt[plta:pltb,plta:pltb]*1e3,w_range,cmap=cm.seismic)
plt.colorbar(label='$\overline{W} \ [mm/s]$', format='%1.3f')
plt.text(50,120,'$W_{max}=$ %1.3f $mm/s$' % (wekmax))
plt.text(50,50,'$W_{min}=$ %1.3f $mm/s$' % (wekmin))
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$U_e$ calculated from u - u(z=300)')
#
plt.suptitle('5 days averaged Ekman pumping $\overline{W_e}$, day %d-%d' % (dstart, dend))
plt.tight_layout(pad=1)
plt.savefig('./figures/W_ekman%d_%d.png' % (dstart, dend))
#============================================================
#Ekman transport Mx
fig;
ax1 = fig.add_subplot(1, 2, 1)
plt.contourf(xc_dom,yc_dom,Mxm[plta:pltb, plta:pltb],100,cmap=cm.seismic)
plt.colorbar(label='$\overline{M_x} \ [m2/s]$', format='%1.3f')
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$U_e = U_t - U_g$ where $U_g$ calculated from geostrophic balance')
#
ax2 = fig.add_subplot(1, 2, 2)
plt.contourf(xc_dom,yc_dom,Mxm_wt[plta:pltb, plta:pltb],100,cmap=cm.seismic)
plt.colorbar(label='$\overline{M_x} \ [m2/s]$', format='%1.3f')
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$U_e$ calculated from u - u(z=300)')
#
plt.suptitle('5 days averaged Ekman transport $\overline{M_x}$, day %d-%d' % (dstart, dend))
plt.tight_layout(pad=1)
plt.savefig('./figures/Mx_day%d_%d.png' % (dstart, dend))

#Ekman transport My
fig;
ax1 = fig.add_subplot(1, 2, 1)
plt.contourf(xc_dom,yc_dom,Mym[plta:pltb, plta:pltb],100,cmap=cm.seismic)
plt.colorbar(label='$\overline{M_x} \ [m2/s]$', format='%1.3f')
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$V_e = V_t - V_g$ where $V_g$ calculated from geostrophic balance')
#
ax2 = fig.add_subplot(1, 2, 2)
plt.contourf(xc_dom,yc_dom,Mym_wt[plta:pltb, plta:pltb],100,cmap=cm.seismic)
plt.colorbar(label='$\overline{M_x} \ [m2/s]$', format='%1.3f')
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$V_e$ calculated from v - v(z=300)')
#
plt.suptitle('5 days averaged Ekman transport $\overline{M_y}$, day %d-%d' % (dstart, dend))
plt.tight_layout(pad=1)
plt.savefig('./figures/My_day%d_%d.png' % (dstart, dend))

# Kinetic Energy 
plt.figure(figsize=(12,6))
plt.plot(time,KEg[:,int(nx/2+20),int(nx/2+20)], label='geostrophic')
plt.plot(time,KEag[:,int(nx/2+20),int(nx/2+20)], label='ageostrophic')
plt.ylabel(r'$KE$')
plt.xlabel("time (hour)")
plt.legend()
plt.title(r'$KE$ at 25km from eddy center, day %d-%d' % (dstart, dend))
plt.savefig('./figures/KE.png')

#
#
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
CS1 = plt.contour(time[hstart1:hend1].squeeze()*2,YC[:,125]*1e-3,Wall[33,:,125,hstart1:hend1]*1e3, levels, colors='0.6')
plt.clabel(CS1, fmt='%2.3f', colors='k', fontsize=10)
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
