# Compute Ekman transport and Ekman pumping from Wenegrat & Thomas 2017

def Me_wenegrat_thomas():
# Compute Ekman transport Wenegrat & Thomas
    dim_factor = np.max(taux)/(rho0*f0)
    Mx_WT = dim_factor*Ro*(Vort_cnd-2*Omega)*np.sin(thetaC)*np.cos(thetaC)/((1+Ro*2*Omega)*(1+Ro*Vort_cnd)-Ro**2*Omega**2)
    My_WT = -dim_factor*(1+Ro*Omega+Ro*2*Omega*np.sin(thetaC)**2+Ro*Vort_cnd*np.cos(thetaC)**2)/((1+Ro*2*Omega)*(1+Ro*Vort_cnd)-Ro**2*Omega**2)
#   In the limit of small Rossby number, simplified formulas are written as
#    Mx_WT = (np.max(taux)/(rho0*f0))*Ro*(Vort_cnd-2*Omega)*np.sin(thetaC)*np.cos(thetaC)
#    My_WT = -(np.max(taux)/(rho0*f0))*(1-Ro*((Vort_cnd-Omega)*np.sin(thetaC)**2 + Omega*np.cos(thetaC)**2))
    return Mx_WT, My_WT
