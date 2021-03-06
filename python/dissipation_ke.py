
"""
Python script to compare dissipation in terms of kinetic energy
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as colors

plt.ion()

def vorticity(U,V,dXC,dYC):
    nx = U[0,:].size
    ny = U[:,0].size
#    nr = U[:,0,0].size
    dyU  = np.zeros((ny,nx));
    dxV  = np.zeros((ny,nx));
    zeta = np.zeros((ny,nx));
    for i in range(1,nx):
        for j in range(1,ny):
            dyU[j,i] = (U[j,i]-U[j-1,i])/dYC[j,i];
            dxV[j,i] = (V[j,i]-V[j,i-1])/dXC[j,i];
            zeta[j,i] = dxV[j,i]-dyU[j,i];
    return zeta

def kinetic_energy(U,V,XC,YC,XG,YG):
    U_c = 0.0*U
    V_c = 0.0*V
    nx = U[0,0,:].size
    ny = U[0,:,0].size
    nr = U[:,0,0].size
# Must interpolate to cell center
    for z in np.arange(0,nr):
        for y in np.arange(0,ny):
            U_c[z,y,:] = np.interp(XC[y,:],XG[y,:],U[z,y,:])
        for x in np.arange(0,nx):
            V_c[z,:,x] = np.interp(YC[:,x],YG[:,x],V[z,:,x])
    KE = U_c**2+V_c**2
    return KE

# select ke domain
pltabc = 80
pltbbc = 250
idepthbc = 21
posextbc = 400
depth_plot_ratio_bc= 0.7# 7

pltabt = 50
pltbbt = 130
idepthbt = 17
posextbt = 500
depth_plot_ratio_bt= 0.9# 7


# Run folders to compare
"""
44G     ./run50a_bc_a25_0
44G     ./run50c_bc_c25_0
86G     ./run51a_bc_a25_ua_Av1e-2
128G    ./run51b_bc_a25_ua_kpp
86G     ./run52a_bc_c25_ua_Av1e-2
66G     ./run52b_bc_c25_ua_kpp
8.6G    ./run60a_bt_a25_0_Av1e-2
8.6G    ./run60c_bt_c25_0_Av1e-2
26G     ./run61a_bt_a25_ua_av1e-2_temp
19G     ./run61b_bt_a25_ua_kpp_sst
26G     ./run62a_bt_c25_ua_av1e-2_sst
19G     ./run62b_bt_c25_ua_kpp_sst
26G     ./run63_bts_a25_ua_Av1e-2
44G     ./run71_bc_a25_ua_Av1e-2
44G     ./run72_bc_c25_ua_Av1e-2
48G     ./run80_bc_a25_kpp_ptracer
62G     ./run81b_bc_a25_kpp_ptracer
48G     ./run81_bc_a25_kpp_ptracer_test_relax
71G     ./run82b_bc_c25_kpp_ptracer

run0 = './run60a_bt_a25_0_Av1e-2/'
run1 = './run61a_bt_a25_ua_av1e-2_temp/'
run2 = './run61b_bt_a25_ua_kpp_sst/'
run3 = './run50a_bc_a25_0/'
run4 = './run51a_bc_a25_ua_Av1e-2/'
run5 = './run51b_bc_a25_ua_kpp/'

run0 = './run60c_bt_c25_0_Av1e-2/'
run1 = './run62a_bt_c25_ua_av1e-2_sst/'
run2 = './run62b_bt_c25_ua_kpp_sst/'
run3 = './run50c_bc_c25_0/'
run4 = './run52a_bc_c25_ua_Av1e-2/'
run5 = './run52b_bc_c25_ua_kpp/'

"""
run0 = './run60a_bt_a25_0_Av1e-2/'
run1 = './run60c_bt_c25_0_Av1e-2/'
run2 = './run61a_bt_a25_ua_av1e-2_temp/'
run3 = './run62a_bt_c25_ua_av1e-2_sst/'
run4 = './run50a_bc_a25_0/'
run5 = './run50c_bc_c25_0/'
run6 = './run51a_bc_a25_ua_Av1e-2/'
run7 = './run52a_bc_c25_ua_Av1e-2/'
run8 = './run51b_bc_a25_ua_kpp/'
run9 = './run52b_bc_c25_ua_kpp/'


run0_label =r'BT_ref acy' #
run1_label =r'BT_ref cy' #
run2_label =r'BT_Av_1e-2 acy' #'Cyclone aut $\nu_h = 10, \nu_v =4e-3 m^2/s$' #
run3_label =r'BT_Av_1e-2 cy' #'Cyclone $\nu_h = 10, \nu_v =1e-5 m^2/s$'
run4_label =r'BC_ref acy' # 'Cyclone $\nu_h = 8, \nu_v =1e-5 m^2/s$' 
run5_label =r'BC_ref cy' #
run6_label =r'BC_Av_1e-2 acy' #'Anticyclone $\nu_h = 8, \nu_v =1e-5 m^2/s$'
run7_label =r'BC_Av_1e-2 cy' #
run8_label =r'BC_KPP acy' #'Anticyclone $\nu_h = 8, \nu_v =1e-5 m^2/s$'
run9_label =r'BC_KPP cy' #

#vort_range = np.linspace(-0.55, 0.55, 101, endpoint=True)
#vorticks = np.linspace(-0.55, 0.55, 11, endpoint=True)
#zero_contour = [0]

XCbt  = mit.rdmds(run0+'XC')
YCbt  = mit.rdmds(run0+'YC')
RCbt  = mit.rdmds(run0+'RC')
RCbc = mit.rdmds(run4+'RC')
XCbc  = mit.rdmds(run4+'XC')
YCbc  = mit.rdmds(run4+'YC')

XGbt  = mit.rdmds(run0+'XG')
YGbt  = mit.rdmds(run0+'YG')
XGbc  = mit.rdmds(run4+'XG')
YGbc  = mit.rdmds(run4+'YG')


nrbt = RCbt.size
nxbt = XCbt[0,:].size
nybt = YCbt[:,0].size

nrbc = RCbc.size
nxbc = XCbc[0,:].size
nybc = YCbc[:,0].size

timestepbc = 120 #timestep input in second
dumpfreqbc = 7200 # in second

timestepbt = 180 #timestep input in second
dumpfreqbt = 10800 # in second

day_s = 0
day_e = 30 #10
day_w = 30 #20
tsbt = dumpfreqbt/timestepbt  # time step = ts*dt (in second); = 7200 = dumpfreq
tsbc = dumpfreqbc/timestepbc 

startfile_bt = day_s*tsbt #or any integer*file_dit
endfile_bt = day_e*86400/timestepbt + 50 #1 day 
endfile_btw = day_w*86400/timestepbt + 50 #1 day 

startfile_bc = day_s*tsbc #or any integer*file_dit
endfile_bc = day_e*86400/timestepbc + 50 #1 day 
endfile_bcw = day_w*86400/timestepbc + 50 #1 day 


itrsbt = np.arange(startfile_bt,endfile_bt, tsbt)
nitbt=itrsbt.size
timebt = itrsbt*timestepbt/3600

itrsbtw = np.arange(startfile_bt,endfile_btw, tsbt)
nitbtw=itrsbtw.size
timebtw = itrsbtw*timestepbt/3600

itrsbc = np.arange(startfile_bc,endfile_bc, tsbc)
nitbc=itrsbc.size
timebc = itrsbc*timestepbc/3600

itrsbcw = np.arange(startfile_bc,endfile_bcw, tsbc)
nitbcw=itrsbcw.size
timebcw = itrsbcw*timestepbc/3600

#icx = int(nx/2)
#icy = int(ny/2)

KE0i = np.zeros((nitbt))
KE1i = np.zeros((nitbt))
KE2i = np.zeros((nitbt))
KE3i = np.zeros((nitbt))
KE4i = np.zeros((nitbc))
KE5i = np.zeros((nitbc))
KE6i = np.zeros((nitbc))
KE7i = np.zeros((nitbc))
KE8i = np.zeros((nitbc))
KE9i = np.zeros((nitbc))

KE0i51 = np.zeros((nitbt))
KE1i51 = np.zeros((nitbt))
KE2i51 = np.zeros((nitbt))
KE3i51 = np.zeros((nitbt))
KE4i51 = np.zeros((nitbc))
KE5i51 = np.zeros((nitbc))
KE6i51 = np.zeros((nitbc))
KE7i51 = np.zeros((nitbc))
KE8i51 = np.zeros((nitbc))
KE9i51 = np.zeros((nitbc))


rho0 = 999.8
f0 = 6.5e-5
f00 = 1e-4

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


#Plot surface vorticity
"""
fig1 = plt.figure(figsize=(4.5,3))
figt = fig1.add_subplot(1, 1, 1)
plt.xlabel("x (km)")
plt.ylabel(r'$\frac{\zeta_0}{|f|}$')
plt.title(r'Surface vorticity section')
###
"""
for it in itrsbt:
    U0 = mit.rdmds(run0+'U',it)
    V0 = mit.rdmds(run0+'V',it)
    KE0 = kinetic_energy(U0,V0,XCbt, YCbt, XGbt, YGbt)
    KE0i[int((it-itrsbt[0])/tsbt)] = np.sum(KE0[:,pltabt:pltbbt,pltabt:pltbbt])
    KE0i51[int((it-itrsbt[0])/tsbt)] = np.sum(KE0[:idepthbt,pltabt:pltbbt,pltabt:pltbbt])
#
    U1 = mit.rdmds(run1+'U',it)
    V1 = mit.rdmds(run1+'V',it)
    KE1 = kinetic_energy(U1,V1,XCbt, YCbt, XGbt, YGbt)
    KE1i[int((it-itrsbt[0])/tsbt)] = np.sum(KE1[:,pltabt:pltbbt,pltabt:pltbbt])
    KE1i51[int((it-itrsbt[0])/tsbt)] = np.sum(KE1[:idepthbt,pltabt:pltbbt,pltabt:pltbbt])
#
    U2 = mit.rdmds(run2+'U',it)
    V2 = mit.rdmds(run2+'V',it)
    KE2 = kinetic_energy(U2,V2,XCbt, YCbt, XGbt, YGbt)
    KE2i[int((it-itrsbt[0])/tsbt)] = np.sum(KE2[:,pltabt:pltbbt,pltabt:pltbbt])
    KE2i51[int((it-itrsbt[0])/tsbt)] = np.sum(KE2[:idepthbt,pltabt:pltbbt,pltabt:pltbbt])

    U3 = mit.rdmds(run3+'U',it)
    V3 = mit.rdmds(run3+'V',it)
    KE3 = kinetic_energy(U3,V3,XCbt, YCbt, XGbt, YGbt)
    KE3i[int((it-itrsbt[0])/tsbt)] = np.sum(KE3[:,pltabt:pltbbt,pltabt:pltbbt])
    KE3i51[int((it-itrsbt[0])/tsbt)] = np.sum(KE3[:idepthbt,pltabt:pltbbt,pltabt:pltbbt])


"""
for it in itrsbtw:
    U1 = mit.rdmds(run1+'U',it)
    V1 = mit.rdmds(run1+'V',it)
#    U2 = mit.rdmds(run2+'U',it)
#    V2 = mit.rdmds(run2+'V',it)
#
    zeta1 = vorticity(U1, V1, dXCbt, dYCbt)
#    zeta2 = vorticity(U2, V2, dXCbt, dYCbt)
#
    svortc1[int((it-itrsbtw[0])/tsbt)]=np.max(abs(zeta1[0,:,:]))
#    svortc2[int((it-itrsbtw[0])/tsbt)]=np.max(abs(zeta2[0,:,:]))
#
#    figt.plot(XC[sec_y,plta:pltb]*1e-3,zeta[0,sec_y,plta:pltb]/f0, label='t=%d hr' %(int(it*timestep/3600)))
#
"""
for it in itrsbc:
    U4 = mit.rdmds(run4+'U',it)
    V4 = mit.rdmds(run4+'V',it)
    KE4 = kinetic_energy(U4,V4,XCbc, YCbc, XGbc, YGbc)
    KE4i[int((it-itrsbc[0])/tsbc)] = np.sum(KE4[:,pltabc:pltbbc,pltabc:pltbbc])
    KE4i[int((it-itrsbc[0])/tsbc)] = np.sum(KE4[:,pltabc:pltbbc,pltabc:pltbbc])

#
    U5 = mit.rdmds(run5+'U',it)
    V5 = mit.rdmds(run5+'V',it)
    KE5 = kinetic_energy(U5,V5,XCbc, YCbc, XGbc, YGbc)
    KE5i[int((it-itrsbc[0])/tsbc)] = np.sum(KE5[:,pltabc:pltbbc,pltabc:pltbbc])
    KE5i51[int((it-itrsbc[0])/tsbc)] = np.sum(KE5[:idepthbc,pltabc:pltbbc,pltabc:pltbbc])

#
    U6 = mit.rdmds(run6+'U',it)
    V6 = mit.rdmds(run6+'V',it)
    KE6 = kinetic_energy(U6,V6,XCbc, YCbc, XGbc, YGbc)
    KE6i[int((it-itrsbc[0])/tsbc)] = np.sum(KE6[:,pltabc:pltbbc,pltabc:pltbbc])
    KE6i51[int((it-itrsbc[0])/tsbc)] = np.sum(KE6[:idepthbc,pltabc:pltbbc,pltabc:pltbbc])

#
    U7 = mit.rdmds(run7+'U',it)
    V7 = mit.rdmds(run7+'V',it)
    KE7 = kinetic_energy(U7,V7,XCbc, YCbc, XGbc, YGbc)
    KE7i[int((it-itrsbc[0])/tsbc)] = np.sum(KE7[:,pltabc:pltbbc,pltabc:pltbbc])
    KE7i51[int((it-itrsbc[0])/tsbc)] = np.sum(KE7[:idepthbc,pltabc:pltbbc,pltabc:pltbbc])

#
    U8 = mit.rdmds(run8+'U',it)
    V8 = mit.rdmds(run8+'V',it)
    KE8 = kinetic_energy(U8,V8,XCbc, YCbc, XGbc, YGbc)
    KE8i[int((it-itrsbc[0])/tsbc)] = np.sum(KE8[:,pltabc:pltbbc,pltabc:pltbbc])
    KE8i51[int((it-itrsbc[0])/tsbc)] = np.sum(KE8[:idepthbc,pltabc:pltbbc,pltabc:pltbbc])

#
    U9 = mit.rdmds(run9+'U',it)
    V9 = mit.rdmds(run9+'V',it)
    KE9 = kinetic_energy(U9,V9,XCbc, YCbc, XGbc, YGbc)
    KE9i[int((it-itrsbc[0])/tsbc)] = np.sum(KE9[:,pltabc:pltbbc,pltabc:pltbbc])
    KE9i51[int((it-itrsbc[0])/tsbc)] = np.sum(KE9[:idepthbc,pltabc:pltbbc,pltabc:pltbbc])

"""
for it in itrsbcw:
    U4 = mit.rdmds(run4+'U',it)
    V4 = mit.rdmds(run4+'V',it)
#    U5 = mit.rdmds(run5+'U',it)
#    V5 = mit.rdmds(run5+'V',it)
#
    zeta4 = vorticity(U4, V4, dXCbc, dYCbc)
#    zeta5 = vorticity(U5, V5, dXCbc, dYCbc)
#
    svortc4[int((it-itrsbcw[0])/tsbc)]=np.max(abs(zeta4[0,:,:]))
#    svortc5[int((it-itrsbcw[0])/tsbc)]=np.max(abs(zeta5[0,:,:]))
#
"""
# surface vort core
fig = plt.figure(figsize=(7,5))
plt.plot(timebt,KE0i, color='k',label=run0_label)
plt.plot(timebt,KE1i, color='k',ls='--', label=run1_label)
plt.plot(timebt,KE2i, color='b',label=run2_label)
plt.plot(timebt,KE3i, color='b',ls='--', label=run3_label)
plt.plot(timebc,KE4i, color='0.6', label=run4_label)
plt.plot(timebc,KE5i, color='0.6',ls='--', label=run5_label)
plt.plot(timebc,KE6i, color='c', label=run6_label)
plt.plot(timebc,KE7i, color='c',ls='--', label=run7_label)
plt.plot(timebc,KE8i, color='r', label=run8_label)
plt.plot(timebc,KE9i, color='r',ls='--', label=run9_label)

plt.xlabel("time (hour)")
plt.ylabel(r'$\int_v KE$')
plt.grid()
plt.legend(loc='best', fancybox=True, framealpha=0.5)
plt.xlim(0,itrsbt[-1]*timestepbt/3600)
#plt.ylim(0.25,0.55)
plt.title(r'Total Kinetic Energy$')
plt.tight_layout(pad=1)
plt.savefig("./figures/dissipation_ke_all.png")

fig = plt.figure(figsize=(7,5))
plt.semilogy(timebt,KE0i, color='k',label=run0_label)
plt.semilogy(timebt,KE1i, color='k',ls='--', label=run1_label)
plt.semilogy(timebt,KE2i, color='b',label=run2_label)
plt.semilogy(timebt,KE3i, color='b',ls='--', label=run3_label)
plt.semilogy(timebc,KE4i, color='0.6', label=run4_label)
plt.semilogy(timebc,KE5i, color='0.6',ls='--', label=run5_label)
plt.semilogy(timebc,KE6i, color='c', label=run6_label)
plt.semilogy(timebc,KE7i, color='c',ls='--', label=run7_label)
plt.semilogy(timebc,KE8i, color='r', label=run8_label)
plt.semilogy(timebc,KE9i, color='r',ls='--', label=run9_label)

plt.xlabel("time (hour)")
plt.ylabel(r'$\int_v KE$')
plt.grid()
plt.legend(loc='best', fancybox=True, framealpha=0.5)
plt.xlim(0,itrsbt[-1]*timestepbt/3600)
#plt.ylim(0.25,0.55)
plt.title(r'Total Kinetic Energy$')
plt.tight_layout(pad=1)
plt.savefig("./figures/dissipation_lnke_all.png")
# ======
fig = plt.figure(figsize=(7,5))
plt.plot(timebt,KE0i51, color='k',label=run0_label)
plt.plot(timebt,KE1i51, color='k',ls='--', label=run1_label)
plt.plot(timebt,KE2i51, color='b',label=run2_label)
plt.plot(timebt,KE3i51, color='b',ls='--', label=run3_label)
plt.plot(timebc,KE4i51, color='0.6', label=run4_label)
plt.plot(timebc,KE5i51, color='0.6',ls='--', label=run5_label)
plt.plot(timebc,KE6i51, color='c', label=run6_label)
plt.plot(timebc,KE7i51, color='c',ls='--', label=run7_label)
plt.plot(timebc,KE8i51, color='r', label=run8_label)
plt.plot(timebc,KE9i51, color='r',ls='--', label=run9_label)

plt.xlabel("time (hour)")
plt.ylabel(r'$\int_v KE$')
plt.grid()
plt.legend(loc='best', fancybox=True, framealpha=0.5)
plt.xlim(0,itrsbt[-1]*timestepbt/3600)
#plt.ylim(0.25,0.55)
plt.title(r'Total Kinetic Energy$')
plt.tight_layout(pad=1)
plt.savefig("./figures/dissipation_ke51.png")

fig = plt.figure(figsize=(7,5))
plt.semilogy(timebt,KE0i51, color='k',label=run0_label)
plt.semilogy(timebt,KE1i51, color='k',ls='--', label=run1_label)
plt.semilogy(timebt,KE2i51, color='b',label=run2_label)
plt.semilogy(timebt,KE3i51, color='b',ls='--', label=run3_label)
plt.semilogy(timebc,KE4i51, color='0.6', label=run4_label)
plt.semilogy(timebc,KE5i51, color='0.6',ls='--', label=run5_label)
plt.semilogy(timebc,KE6i51, color='c', label=run6_label)
plt.semilogy(timebc,KE7i51, color='c',ls='--', label=run7_label)
plt.semilogy(timebc,KE8i51, color='r', label=run8_label)
plt.semilogy(timebc,KE9i51, color='r',ls='--', label=run9_label)

plt.xlabel("time (hour)")
plt.ylabel(r'$\int_v KE$')
plt.grid()
plt.legend(loc='best', fancybox=True, framealpha=0.5)
plt.xlim(0,itrsbt[-1]*timestepbt/3600)
#plt.ylim(0.25,0.55)
plt.title(r'Total Kinetic Energy$')
plt.tight_layout(pad=1)
plt.savefig("./figures/dissipation_lnke51.png")

"""
itrsd = [60,480,960,1440,1920,2400]
fig = plt.figure(figsize=(12,6))
ax1 = plt.subplot(131)
plt.xlabel(r'$\zeta/f_0$')
plt.ylabel("Depth [m]")
plt.title('Vorticity at initial eddy core, case Av = 1e-2')
plt.grid(True)
ax2 = plt.subplot(132, sharex=ax1, sharey=ax1)
plt.xlabel(r'$\zeta/f_0$')
#plt.ylabel("Depth [m]")
plt.title('Vorticity at initial eddy core, case with KPP')
plt.grid(True)
ax3 = plt.subplot(133, sharex=ax1, sharey=ax1)
plt.xlabel(r'$\zeta/f_0$')
#plt.ylabel("Depth [m]")
plt.title('Vorticity at initial eddy core, case Av=4e-3')
plt.grid(True)


for it in itrsd:
    viscKPP3 = mit.rdmds(run3+'KPPviscAz',it)
    viscKPP5 = mit.rdmds(run5+'KPPviscAz',it)
    U = mit.rdmds(run1+'U',it)
    V = mit.rdmds(run1+'V',it)
    zetaz = vorticity(U, V, dXC, dYC)
    Uk = mit.rdmds(run3+'U',it)
    Vk = mit.rdmds(run3+'V',it)
    zetak = vorticity(Uk, Vk, dXC, dYC)
    U100 = mit.rdmds(run4+'U',it)
    V100 = mit.rdmds(run4+'V',it)
    zeta100 = vorticity(U100, V100, dXC, dYC)

###########################################################    
    plt.figure (figsize = (4.5, 6))
    plt.plot(np.tile(1e-5,(RC.size)), RC.squeeze(), label=run0_label)
    plt.plot(np.tile(4e-3,(RC.size)), RC.squeeze(), label=run1_label)
    plt.plot(np.tile(1e-2,(RC.size)), RC.squeeze(), label=run2_label)
    plt.plot(viscKPP3[:,90,90], RC.squeeze(),ls='--', label=run3_label)
    plt.plot(np.tile(4e-3,(RC.size)), RC.squeeze(), label=run4_label)
    plt.plot(viscKPP5[:,90,90], RC200.squeeze(),ls='--', label=run5_label)
#    plt.plot(np.tile(1e-5,(RC.size)), RC.squeeze())
#    plt.plot(viscKPP0[:,90,90], RC.squeeze())
    plt.xlabel(r'$\nu_z [m^2/s]$')
    plt.ylabel("Depth [m]")
    plt.title("viscosity vertical")
    plt.tight_layout(pad=1)
    plt.legend(loc='best')
    plt.savefig('./figures/Compare_Az_%d.png' %(it))
#################################################################
#
    ax1.plot(zetaz[:,90,90]/f0, RC.squeeze(), label='timestep %d hr' % (it*timestep/3600))
    ax2.plot(zetak[:,90,90]/f0, RC.squeeze(), label='timestep %d hr' % (it*timestep/3600))
    ax3.plot(zeta100[:,90,90]/f0, RC.squeeze()/3, label='timestep %d hr' % (it*timestep/3600))

ax3.legend(loc='lower right')
plt.tight_layout(pad=1)
plt.savefig('./figures/Vorticityprofilekppvsno_100.png')
"""
############

"""
# surface vort core difference
fig = plt.figure(figsize=(7,4))
plt.plot(time,-(svortc1-svortc0)/f0, label='anticyclone')
plt.plot(time,(svortc2+svortc0)/f0, label='cyclone')
plt.xlabel("time (hour)")
plt.ylabel(r'$\frac{|\zeta_c|-|\zeta_c^{ref}|}{|f|}$', fontsize=14)
plt.legend()
plt.xlim(0,itrs[-1]*timestep/3600)
plt.title(r'Wind modified core vorticity, $\frac{|\zeta_c|-|\zeta_c^{ref}|}{|f|}$', fontsize=14)
plt.tight_layout(pad=1)
plt.savefig("./figures/modcore.png")

"""
