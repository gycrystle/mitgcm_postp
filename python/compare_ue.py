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
pltb = 150
idepth = 99

timestep = 180 #timestep input in second
dumpfreq = 10800 # in second

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
#itrs=np.arange(nstart,nend,1)
itrs=[-1,0,1,5,10,20,50,70,90,120]#,480,960,1440,1920,2400]
#nit=itrs.size
nit=4

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

Ug300_all=np.zeros((nit,nr,ny,nx))
Vg300_all=np.zeros((nit,nr,ny,nx))

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
W_sternall = np.zeros((int(nit/8+1),ny,nx));

termx2 = np.zeros((ny,nx));
termy2 = np.zeros((ny,nx));
#
dytermx2 = np.zeros((ny,nx));
dxtermy2 = np.zeros((ny,nx));
W_stern2 = np.zeros((ny,nx));
W_sternall2 = np.zeros((int(nit/8+1),ny,nx));


dstart = day_s #5
dend = day_e # (itrs[-1]+1)*ts/720  #10

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
#    P = mit.rdmds('PH',((it+1)*ts))
    U_thta = V*np.cos(thetaC) - U*np.sin(thetaC)
    U_r = U*np.cos(thetaC)+V*np.sin(thetaC)
#
    """
    #dxP = np.diff(P, axis=2)/np.tile(dXC[:,1:],[nr,1,1])
    dxP = np.diff(P, axis=2)/dX_Ug
    Vg = dxP/f0
    dxP300 = np.diff(P[idepth,:,:], axis=1)/np.diff(XC, axis=1)
    Vg300 = dxP300/f0
#
    #dyP = np.diff(P, axis=1)/np.tile(dYC[1:,:],[nr,1,1])
    dyP = np.diff(P, axis=1)/dY_Vg
    Ug = -dyP/f0
    dyP300 = np.diff(P[idepth,:,:], axis=0)/np.diff(YC, axis=0)
    Ug300 = -dyP300/f0
#
# compute U ageostrophic from U bottom
    u_wt = U - U[idepth,:,:]
    v_wt = V - V[idepth,:,:] 
#
#Interpolate Ug in cell center, losing first grids(bottom and left)
    for z in np.arange(0,nr):
        for y in np.arange(0,ny):
            Vg_c[z,y,:] = np.interp(XC[y,:],XUg[y,:-1],Vg[z,y,:])
            Vg300_c[y,:] = np.interp(XC[y,:],XUg[y,:-1],Vg300[y,:])
            Ut_c[z,y,:] = np.interp(XC[y,:],XG[y,:],U[z,y,:])
        for x in np.arange(0,nx):
            Ug_c[z,:,x] = np.interp(YC[:,x],YUg[:-1,x],Ug[z,:,x])
            Ug300_c[:,x] = np.interp(YC[:,x],YUg[:-1,x],Ug300[:,x])
            Vt_c[z,:,x] = np.interp(YC[:,x],YG[:,x],V[z,:,x]) #verify if V points = YUg
#
# transform u_geostrophic to cylindrical to compute cyclogeo vel
    U_theta_g = Vg_c*np.cos(thetaC)-Ug_c*np.sin(thetaC)
    U_r_g = Ug_c*np.cos(thetaC) + Vg_c*np.sin(thetaC)
# Tangential vel component of gradient wind, cell centered
    U_theta_cygeo = (f0*radC/2)-np.sqrt(f0**2*radC**2/4 - f0*radC*U_theta_g)
# transform U_cyclogeo to cartesian
    U_cygeo = - U_theta_cygeo*np.sin(thetaC) # Ur_cygeo*np.cos(thetaC)
    V_cygeo = U_theta_cygeo*np.cos(thetaC) # Ur_cygeo*np.sin(thetaC)
#
    Vg_all[(it-itrs[0]),:,:,:]= Vg_c#np.tile(Vg300_c,(nr,1,1)) #cell centered Vg
    Ug_all[(it-itrs[0]),:,:,:]= Ug_c#np.tile(Ug300_c,(nr,1,1)) #cell centered Ug
    Vg300_all[(it-itrs[0]),:,:,:]= np.tile(Vg300_c,(nr,1,1)) #cell centered Vg
    Ug300_all[(it-itrs[0]),:,:,:]= np.tile(Ug300_c,(nr,1,1)) #cell centered Ug
    Ut_all[(it-itrs[0]),:,:,:]= Ut_c; #cell centered u_total
    Vt_all[(it-itrs[0]),:,:,:]= Vt_c; #cell centered v total
    U_cg_all[(it-itrs[0]),:,:,:]= U_cygeo #np.tile(U_cygeo, (nr,1,1));
    V_cg_all[(it-itrs[0]),:,:,:]= V_cygeo #np.tile(V_cygeo, (nr,1,1));
#
    u_wtall[(it-itrs[0]),:,:,:]= u_wt; #u-point u_wt
    v_wtall[(it-itrs[0]),:,:,:]= v_wt; #v-point v_wt
    """
    plt.figure()
    plt.contourf(xc_dom,yc_dom,U_r[99,plta:pltb,plta:pltb],100,cmap=cm.seismic)#
    plt.colorbar()
#
#======= time loop end

################## Compute ageostrophic flow ##############

#1 from geostrophic balance
Ue_g = Ut_all - Ug_all
Ve_g = Ut_all - Vg_all

#2 from geostrophic balance at z=300
Ue_g300 = Ut_all - Ug300_all
Ve_g300 = Ut_all - Vg300_all

#3 from cyclogeostrophic balance
Ue_cyg = Ut_all - U_cg_all
Ve_cyg = Vt_all - V_cg_all


#============ PLOT ==========================
# Various U at t=0
w_range = np.linspace(-0.4, 0.4, 101, endpoint=True)
#
levels = np.concatenate((np.linspace(-0.5,0,10,endpoint=False),np.linspace(0.05,0.5,10,endpoint=True)),axis=0)
#
time = (itrs)*dumpfreq/3600
#

######################################################################################
# Ekman tranport plot
we_range = np.linspace(-0.4, 0.4, 101, endpoint=True)
fig = plt.figure(figsize=(14,4))
ax1 = plt.subplot(131)
plt.contourf(xc_dom,yc_dom,Ut_all[0,99,plta:pltb,plta:pltb],we_range,cmap=cm.seismic)#
cb=plt.colorbar(format='%1.3f')
cb.set_label(r'$U \, (m/s)$', labelpad=-40, y=1.1, rotation=0)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
ax1.set_aspect(1)
plt.text(480,480,'$U_{extrema}=$ [%1.3f, %1.3f] $m/s$' % (np.min(Ut_all[0,99,plta:pltb,plta:pltb]),np.max(Ut_all[0,99,plta:pltb,plta:pltb])), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title(r'$U (z=-300)$')
#plt.title(r'$\overline{W_{e}}=\overline{\partial_x \int \, U_eg + \partial_y \int \, V_eg}$ day %d-%d' % (dstart, dend),fontsize=11)
#
ax2 = plt.subplot(132, sharex=ax1)
ax2.set_aspect(1)
plt.contourf(xc_dom,yc_dom,Ug300_all[0,99,plta:pltb,plta:pltb],we_range,cmap=cm.seismic)#
cb=plt.colorbar(format='%1.3f')
cb.set_label(r'$U \, (m/s)$', labelpad=-40, y=1.1, rotation=0)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
plt.text(480,480,'$U_{extrema}=$ [%1.3f, %1.3f] $m/s$' % (np.min(Ug300_all[0,99,plta:pltb,plta:pltb]),np.max(Ug300_all[0,99,plta:pltb,plta:pltb])), fontsize=10)
#plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(We_g300)*1e3,np.max(We_g300)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertically averaged W at $z[0, %d]m$, timestep %d hr' % (RC[idepth],int(it*timestep/3600)),fontsize=11)
#plt.title(r'$\overline{W_{e}}=\overline{\partial_x \int \, U_eg300 + \partial_y \int \, V_eg300}$ day %d-%d' % (dstart, dend),fontsize=11)
plt.title(r'$U_{geostrphic} (z=-300)$')
#
ax3 = plt.subplot(133, sharex=ax1)
ax3.set_aspect(1)
plt.contourf(xc_dom,yc_dom,U_cg_all[0,99,plta:pltb,plta:pltb],we_range,cmap=cm.seismic)#
cb=plt.colorbar(format='%1.3f')
cb.set_label(r'$U \, (m/s)$', labelpad=-40, y=1.1, rotation=0)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
#plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(We_wt)*1e3,np.max(We_wt)*1e3), fontsize=10)
plt.text(480,480,'$U_{extrema}=$ [%1.3f, %1.3f] $m/s$' % (np.min(U_cg_all[0,99,plta:pltb,plta:pltb]),np.max(U_cg_all[0,99,plta:pltb,plta:pltb])), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#plt.title(r'$\overline{W_{e}}=\overline{\partial_x \int \, U_e300 + \partial_y \int \, V_e300}$ day %d-%d' % (dstart, dend),fontsize=11)
plt.title(r'$U_{cyclogeo}(z=-300)$')
#
"""
ax4 = plt.subplot(234, sharex=ax1)
ax4.set_aspect(1)
plt.contourf(xc_dom,yc_dom,u_wtall[0,1,plta:pltb,plta:pltb],we_range,cmap=cm.seismic)#
cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
#plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(We_cyg)*1e3,np.max(We_cyg)*1e3), fontsize=10)
plt.text(480,480,'$U_{extrema}=$ [%1.3f, %1.3f] $m/s$' % (np.min(u_wtall),np.max(u_wtall)), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#plt.title(r'$\overline{W_{e}}=\overline{\partial_x \int \, U_{ecg} + \partial_y \int \, V_{ecg}}$ day %d-%d' % (dstart, dend),fontsize=11)

ax5 = plt.subplot(235, sharex=ax1)
ax5.set_aspect(1)
plt.contourf(xc_dom,yc_dom,Ue_g300[0,1,plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
#CS1 = plt.contour(YC[plta:pltb,int(si_x/2)]*1e-3,RC.squeeze(),W[:,plta:pltb,int(si_x/2)]*1e3, levels, colors='0.6')
#plt.clabel(CS1, fmt='%2.2f', colors='k', fontsize=10)
plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(We_stern)*1e3,np.max(We_stern)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#    plt.title('Vertical velocity $W_{num}$ at $x=%dkm$, timestep %d hr' % (XC[int(si_x/2),int(si_x/2)]*1e-3,int(it*timestep/3600)), fontsize=11)
    #
#plt.title(r'$\overline{W_{stern}}=\overline{\frac{1}{\rho_0} \nabla \times \left[\frac{\tau}{f+\zeta}\right]}$', fontsize=14)
###
ax6 = plt.subplot(236, sharex=ax1)
ax6.set_aspect(1)
plt.contourf(xc_dom,yc_dom,Ue_cyg[0,1,plta:pltb,plta:pltb]*1e3,we_range,cmap=cm.seismic)#
cb=plt.colorbar(label="W (mm/s)", format='%1.3f')
cb.set_label(r'$W \, (mm/s)$', labelpad=-40, y=1.1, rotation=0)
#CS2 = plt.contour(xc_dom,yc_dom,Wmint[plta:pltb,plta:pltb]*1e3, levels, colors='0.6')
#plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=10)
#plt.text(480,480,'$W_{extrema}=$ [%1.3f, %1.3f] $mm/s$' % (np.min(We_stern2)*1e3,np.max(We_stern2)*1e3), fontsize=10)
plt.xlabel("x (km)")
plt.ylabel("y (km)")
#plt.title(r'$\overline{W_{stern*}}=\overline{\frac{1}{\rho_0} \nabla \times \left[\frac{\tau}{f+2\zeta}\right]}$', fontsize=14)
"""
plt.tight_layout (pad = 1)
plt.savefig('./figures/U_comp.png')
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
