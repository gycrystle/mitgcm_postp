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
