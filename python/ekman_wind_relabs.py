# Python
# To compare the the ekman pumping cross section for the simulation 
# With relative wind and absolute wind

import numpy as np
import MITgcmutils as mit
from MITgcmutils import rdmds
import matplotlib.pyplot as plt
from matplotlib import cm

plt.ion()

rho0 = 1028.0
XC = mit.rdmds('./run_relativewind/XC')
dXC = mit.rdmds('./run_relativewind/DXC')
YC = mit.rdmds('./run_relativewind/YC')
dYC = mit.rdmds('./run_relativewind/DYC')
RC = mit.rdmds('./run_relativewind/RC')

dyUr = 0.0*np.ones((250,250));
dxVr = 0.0*np.ones((250,250));
dyUa = 0.0*np.ones((250,250));
dxVa = 0.0*np.ones((250,250));


zetar = 0.0*np.ones((250,250));
zetar_intx = 0.0*np.ones((250,250));
zetar_inty = 0.0*np.ones((250,250));
zetaa = 0.0*np.ones((250,250));
zetaa_intx = 0.0*np.ones((250,250));
zetaa_inty = 0.0*np.ones((250,250));
#
termxr =0.0*np.ones((250,250));
termyr =0.0*np.ones((250,250));
termxa =0.0*np.ones((250,250));
termya =0.0*np.ones((250,250));
#
dytermxr = 0.0*np.ones((250,250));
dxtermyr = 0.0*np.ones((250,250));
W_totr = 0.0*np.ones((250,250));
#
dytermxa = 0.0*np.ones((250,250));
dxtermya = 0.0*np.ones((250,250));
W_tota = 0.0*np.ones((250,250));

w_range = np.linspace(-0.5, 0.5, 101, endpoint=True)

itrs = [720, 1440, 2160, 2880, 3600] #, 2160, 2880, 3600, 4320, 5040, 5760, 6480, 7200]


for it in itrs:
    Ur = mit.rdmds('./run_relativewind/U',it)
    Vr = mit.rdmds('./run_relativewind/V',it)
    Ua = mit.rdmds('./run_wdiag_5/U',it)
    Va = mit.rdmds('./run_wdiag_5/V',it)
#
    tauxr = mit.rdmds('./run_relativewind/diagTAUX',it) 
    tauyr = mit.rdmds('./run_relativewind/diagTAUY',it)
#
    tauxa = mit.rdmds('./run_wdiag_5/diagTAUX',it)                                                                                  
    tauya = mit.rdmds('./run_wdiag_5/diagTAUY',it)

    for i in range(1,250):
        for j in range(1,250):
            dyUr[j,i] = (Ur[1,j,i]-Ur[1,j-1,i])/dYC[j,i];
            dxVr[j,i] = (Vr[1,j,i]-Vr[1,j,i-1])/dXC[j,i];
            zetar[j,i] = dxVr[j,i]-dyUr[j,i];
#
            dyUa[j,i] = (Ua[1,j,i]-Ua[1,j-1,i])/dYC[j,i];
            dxVa[j,i] = (Va[1,j,i]-Va[1,j,i-1])/dXC[j,i];
            zetaa[j,i] = dxVa[j,i]-dyUa[j,i];
#
    for i in range(1,249):
        for j in range(1,249):
            zetar_intx[j,i] =(zetar[j,i]+zetar[j+1,i])/2
            zetar_inty[j,i] =(zetar[j,i]+zetar[j,i+1])/2
#
            termxr[j,i] =tauxr[j,i]/(1e-4+zetar_intx[j,i])
            termyr[j,i] =tauyr[j,i]/(1e-4+zetar_inty[j,i])
#
            zetaa_intx[j,i] =(zetaa[j,i]+zetaa[j+1,i])/2
            zetaa_inty[j,i] =(zetaa[j,i]+zetaa[j,i+1])/2
#
            termxa[j,i] =tauxa[j,i]/(1e-4+zetaa_intx[j,i])
            termya[j,i] =tauya[j,i]/(1e-4+zetaa_inty[j,i])
#
    for i in range(1,248):
        for j in range(1,248):
            dytermxr[(j+1),i] = (termxr[j+1,i]-termxr[j,i])/dYC[j+1,i];
            dxtermyr[j,(i+1)] = (termyr[j,i+1]-termyr[j,i])/dXC[j,i+1];
            W_totr[(j+1),(i+1)] = (dxtermyr[j+1,i+1]-dytermxr[j+1,i+1])/rho0;
#
            dytermxa[(j+1),i] = (termxa[j+1,i]-termxa[j,i])/dYC[j+1,i];
            dxtermya[j,(i+1)] = (termya[j,i+1]-termya[j,i])/dXC[j,i+1];
            W_tota[(j+1),(i+1)] = (dxtermya[j+1,i+1]-dytermxa[j+1,i+1])/rho0;


    plt.figure()

    plt.plot(YC[:,125]*1e-3,W_totr[:,125]*1e3,'b', label='Relative wind')
    plt.plot(YC[:,125]*1e-3,W_tota[:,125]*1e3,'r', label='Absolute wind') 
    plt.xlabel("y (km)")
    plt.ylabel("Wtot (mm/s)")
    plt.legend()
    plt.title("Ekman pumping at x=150km, it = "+ str(it))
    
    plt.savefig("ekman_pumping" + str(it) + ".png")

   
