import MITgcmutils as mit
from MITgcmutils import rdmds
import def_radius
import numpy as np

folder = './run21_00btp/'
RC = mit.rdmds(folder+'RC')

#H = np.array([500, 3000]) # Layer thickness (m)
H = np.diff(RC)

temp_ref_argo=sio.loadmat('temp_ref_argo.mat')
temp_ref = temp_ref_argo['temp']
temp_ref_art=sio.loadmat('art_temp_ref.mat')
temp_ref1 = temp_ref_art['mPtempout']

gp = np.array([1e-2])     # reduced gravity (m/s^2) -- density jumps between layers
f = 1e-4                  # coriolis parameter (1/s)

rd = def_radius.cal_rad(H,gp,f)

print("first deformation radius: {0:.1f} km".format(rd[1]*1e-3))
