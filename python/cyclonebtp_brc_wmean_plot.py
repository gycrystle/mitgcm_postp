"""
Plot the wmean on the surface and the cross section
"""

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

plt.ion()

dir0 = './run62a_bt_c25_ua_av1e-2_sst/'
dir1 = './run62b_bt_c25_ua_kpp_sst/'
dir2 = './run52a_bc_c25_ua_Av1e-2/'
dir3 = './run52b_bc_c25_ua_kpp/'

"""
dir0 = './run61a_bt_a25_ua_av1e-2_temp/'
dir1 = './run61b_bt_a25_ua_kpp_sst/'
dir2 = './run51a_bc_a25_ua_Av1e-2/'
dir3 = './run51b_bc_a25_ua_kpp/'
"""
#v_range = 100
v_range_max = 6.5
v_range = np.linspace(-v_range_max, v_range_max, 101, endpoint=True) #.02 or .12 
v_ticks = np.linspace(-v_range_max, v_range_max, 11, endpoint=True)
#levels = np.linspace(-0.02, 0.02,6, endpoint=True)
levels = np.concatenate((np.linspace(-10.0,0,10,endpoint=False),np.linspace(1.0,10.0,10,endpoint=True)),axis=0)
conv_factor = 86400
conv_unit = 'm/day'

day_s = 10
day_e = 15

XC_btp = mit.rdmds(dir0+'XC')
YC_btp = mit.rdmds(dir0+'YC')
RC_btp = mit.rdmds(dir0+'RC')
nx_btp = XC_btp.shape[1]
ny_btp = YC_btp.shape[0]
nr_btp = RC_btp.size
icx_btp = int(nx_btp/2)
icy_btp = int(ny_btp/2)
#idepth = 17 #17 40

# =======  Barotropic ============
#
timestep = 180 #timestep input in second
dumpfreq = 10800 # file every 3 hours for dumpfreq 10800
file_dit = dumpfreq/timestep

startfile = day_s*file_dit #or any integer*file_dit
endfile = day_e*86400/timestep + 50 #1 day 
itrs = np.arange(startfile,endfile, file_dit)
nit_btp = itrs.size

W5d1=np.zeros((nit_btp,nr_btp,ny_btp,nx_btp))
W5d2=np.zeros((nit_btp,nr_btp,ny_btp,nx_btp))

# time loop
for it in itrs:
    W1 = mit.rdmds(dir0+'W',it)
    W5d1[int((it-itrs[0])/file_dit), :,:,:]= W1;
    W2 = mit.rdmds(dir1+'W',it)
    W5d2[int((it-itrs[0])/file_dit), :,:,:]= W2;

W_btp1 = np.mean(W5d1, axis=0)
W_btp2 = np.mean(W5d2, axis=0)
#wmin = np.min(W_ta)*conv_factor
#wmax = np.max(W_ta)*conv_factor
itrs=[]

# =======  Baroclinic ============
#
timestep = 120 #timestep input in second
dumpfreq = 7200 # file every 3 hours for dumpfreq 10800
file_dit = dumpfreq/timestep

startfile = day_s*file_dit #or any integer*file_dit
endfile = day_e*86400/timestep + 50 #1 day 
itrs = np.arange(startfile,endfile, file_dit)
nit_brc = itrs.size

XC_brc = mit.rdmds(dir2+'XC')
YC_brc = mit.rdmds(dir2+'YC')
RC_brc = mit.rdmds(dir2+'RC')
nx_brc = XC_brc.shape[1]
ny_brc = YC_brc.shape[0]
nr_brc = RC_brc.size
icx_brc = int(nx_brc/2)
icy_brc = int(ny_brc/2)

W5d1=np.zeros((nit_brc,nr_brc,ny_brc,nx_brc))
W5d2=np.zeros((nit_brc,nr_brc,ny_brc,nx_brc))

# time loop
for it in itrs:
    W1 = mit.rdmds(dir2+'W',it)
    W5d1[int((it-itrs[0])/file_dit), :,:,:]= W1;
    W2 = mit.rdmds(dir3+'W',it)
    W5d2[int((it-itrs[0])/file_dit), :,:,:]= W2;

W_brc1 = np.mean(W5d1, axis=0)
W_brc2 = np.mean(W5d2, axis=0)


y_btp = np.linspace(0.0,300.0, 120, endpoint=True)
y_brc = np.linspace(0.0,300.0, 250, endpoint=True)

plta = 40
pltb = 290
pltc = 30
pltd = 150

plt.figure(figsize=(5,8))
plt.plot(W_btp1[17,pltc:pltd,90]*conv_factor, y_btp, color='b', label='BTc_1')
plt.plot(W_btp2[17,pltc:pltd,90]*conv_factor, y_btp, color='b', ls='dashed', label='BTc_2')
plt.plot(W_brc1[25,plta:pltb,165]*conv_factor, y_brc, color='r', label='BCc_1')
plt.plot(W_brc2[25,plta:pltb,165]*conv_factor, y_brc, color='r', ls='dashed', label='BCc_2')

plt.legend(loc='best')
plt.grid()
plt.xlabel(r'$\overline{W}$ [m/day]')
plt.ylabel('y (km)')
plt.title(r'Time averaged vertical velocity $\overline{W}$ at z=~51m')
plt.savefig("./figures/cyclone_btp_brc.png")

# vertical section 
#
plt.figure(figsize=(7,6))
ax1 = plt.subplot(121)
plt.xlabel(r'$\overline{W} \ [m/day]$')
plt.ylabel("Depth [m]")
#plt.title('Inside eddy, %d km from center' % (abs(XC[94,94]-XC[90,90])*1e-3))
plt.title('Inside eddy, %d km from center' % (abs(XC_brc[173,173]-XC_brc[165,165])*1e-3))
plt.grid(True)
ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
plt.xlabel(r'$\overline{W} \ [m/day]$')
plt.grid(True)
plt.ylim(-1500,0)
#plt.ylabel("Depth [m]")
#plt.title('Inside eddy, %d km from center' % (abs(XC[104,104]-XC[90,90])*1e-3))
plt.title('Inside eddy, %d km from center' % (abs(XC_brc[194,194]-XC_brc[165,165])*1e-3))

ax1.plot(W_btp1[:,94,90]*conv_factor, RC_btp[:].squeeze(), color='b')#, label='timestep %d hr' % (it*timestep/3600))
ax1.plot(W_btp2[:,94,90]*conv_factor, RC_btp[:].squeeze(), color='b', ls='dashed')#, label='timestep %d hr' % (it*timestep/3600))
ax1.plot(W_brc1[:,173,165]*conv_factor, RC_brc[:].squeeze(), color='r')#, label='timestep %d hr' % (it*timestep/3600))
ax1.plot(W_brc2[:,173,165]*conv_factor, RC_brc[:].squeeze(), color='r', ls='dashed')#, label='timestep %d hr' % (it*timestep/3600))

ax2.plot(W_btp1[:,104,90]*conv_factor, RC_btp[:].squeeze(), color='b', label='BTc_1')
ax2.plot(W_btp2[:,104,90]*conv_factor, RC_btp[:].squeeze(), color='b', ls='dashed', label='BTc_2')
ax2.plot(W_brc1[:,194,165]*conv_factor, RC_brc[:].squeeze(), color='r', label='BCc_1')
ax2.plot(W_brc2[:,194,165]*conv_factor, RC_brc[:].squeeze(), color='r', ls='dashed', label='BCc_2')

plt.tight_layout(pad=1)
ax2.legend(loc='lower right')
plt.savefig('./figures/Wmean_profile_cyclone.png')

plt.figure(figsize=(7,6))
ax1 = plt.subplot(121)
plt.xlabel(r'$\overline{W} \ [m/day]$')
plt.ylabel("Depth [m]")
#plt.title('Inside eddy, %d km from center' % (abs(XC[80,80]-XC[90,90])*1e-3c)
plt.title('Inside eddy, %d km from center' % (abs(XC_brc[145,145]-XC_brc[165,165])*1e-3))
plt.grid(True)
ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
plt.xlabel(r'$\overline{W} \ [m/day]$')
plt.grid(True)
plt.ylim(-1500,0)

plt.title('Inside eddy, %d km from center' % (abs(XC_brc[194,194]-XC_brc[165,165])*1e-3))

ax1.plot(W_btp1[:,80,90]*conv_factor, RC_btp[:].squeeze(), color='b')
ax1.plot(W_btp2[:,80,90]*conv_factor, RC_btp[:].squeeze(), color='b', ls='dashed') #, label='timestep %d hr' % (it*timestep/3600))
ax1.plot(W_brc1[:,145,165]*conv_factor, RC_brc[:].squeeze(), color='r', ls='dashed')#, label='timestep %d hr' % (it*timestep/3600))
ax1.plot(W_brc1[:,145,165]*conv_factor, RC_brc[:].squeeze(), color='r')#, label='timestep %d hr' % (it*timestep/3600))

ax2.plot(W_btp1[:,70,90]*conv_factor, RC_btp[:].squeeze(), color='b', label='BTc_1')
ax2.plot(W_btp2[:,70,90]*conv_factor, RC_btp[:].squeeze(), color='b', ls='dashed', label='BTc_2')
ax2.plot(W_brc1[:,124,165]*conv_factor, RC_brc[:].squeeze(), color='r', label='BCc_1')
ax2.plot(W_brc2[:,124,165]*conv_factor, RC_brc[:].squeeze(), color='r', ls='dashed', label='BCc_2')

plt.tight_layout(pad=1)	
ax2.legend(loc='lower right')
plt.savefig('./figures/Wmean_profiledown_cyclone.png')


