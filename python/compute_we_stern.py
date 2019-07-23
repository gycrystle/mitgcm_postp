# Compute Ekman pumping with Stern formula
#import numpy as np
def we_stern(vort,taux,tauy,dXC,dYC,f0,rho0):
    import numpy as np
    snx = dXC[0,:].size
    sny = dYC[:,0].size
    szeta_intx = np.zeros((sny,snx))
    szeta_inty = np.zeros((sny,snx))
    stermx = np.zeros((sny,snx))
    stermy = np.zeros((sny,snx))
    sdytermx = np.zeros((sny,snx))
    sdxtermy = np.zeros((sny,snx))
    sW_stern = np.zeros((sny,snx))
    for i in range(1,snx-1):
        for j in range(1,sny-1):
            szeta_intx[j,i] =(vort[j,i]+vort[j+1,i])/2 #vorticity @ u point
            szeta_inty[j,i] =(vort[j,i]+vort[j,i+1])/2
            stermx[j,i] =taux[j,i]/(f0+szeta_intx[j,i])
            stermy[j,i] =tauy[j,i]/(f0+szeta_inty[j,i])
    for i in range(1,snx-2):
        for j in range(1,sny-2):
            sdytermx[(j+1),i] = (stermx[j+1,i]-stermx[j,i])/dYC[j+1,i]
            sdxtermy[j,(i+1)] = (stermy[j,i+1]-stermy[j,i])/dXC[j,i+1]
            sW_stern[(j+1),(i+1)] = (sdxtermy[j+1,i+1]-sdytermx[j+1,i+1])/rho0
    return sW_stern
