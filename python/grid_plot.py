import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.patches as patches
#from matplotlib import animation
import MITgcmutils as mit
from MITgcmutils import rdmds

XC = mit.rdmds('XC')
YC = mit.rdmds('YC')

plt.ion()

dx=2.5
dxs = [dx*16,dx*8,dx*4,dx*2,dx,dx*2,dx*4,dx*8,dx*16]
si_xs = [3,3,3,3,10,3,3,3,3]
si_xssum = np.cumsum(si_xs)
#
xx1=np.linspace(0,si_xssum[-1],si_xssum[-1]+1)
xx1[0]=0.0
for ix in np.arange(1,si_xssum[0]+1):
  xx1[ix]=xx1[ix-1]+dxs[0]
for ii in np.arange(1,si_xssum.size):
  for ix in np.arange(si_xssum[ii-1]+1,si_xssum[ii]+1):
    xx1[ix]=xx1[ix-1]+dxs[ii]
#
yy1 = xx1 # ok for square domain
#
dx1 = np.diff(xx1)
dy1 = np.diff(yy1)
#
xx = xx1[0:-1] + 0.5*dx1
yy = yy1[0:-1] + 0.5*dy1
si_x = xx.size
si_y = yy.size

fig = plt.figure(figsize=(9,8))
plt.axis('equal')
#plt.grid()
ax = fig.add_subplot(111)
"""
for i in np.arange(0,si_x):
    plt.plot(np.tile(xx1[i],(si_x+1,1,1)).squeeze(),yy1[:],color='0.7')
    plt.plot(xx1[:],np.tile(yy1[i],(si_x+1,1,1)).squeeze(),color='0.5')

"""
for i in np.arange(0,180):
    plt.plot(XC[i,:]*1e-3,YC[i,:]*1e-3,color='0.7')
    plt.plot(XC[:,i]*1e-3,YC[:,i]*1e-3,color='0.5')


plt.ylim(-50,1210)
plt.xlim(-50,1300)
ax.plot([XC[30,30]*1e-3, XC[150,150]*1e-3, XC[150,150]*1e-3, XC[30,30]*1e-3, XC[30,30]*1e-3], [YC[30,30]*1e-3,YC[30,30]*1e-3, YC[150,150]*1e-3, YC[150,150]*1e-3, YC[30,30]*1e-3], color="blue", linestyle="--")
#ax.plot([xx1[12], xx1[22], xx1[22], xx1[12], xx1[12]], [yy1[12], yy1[12], yy1[22], yy1[22], yy1[12]], color="blue", linestyle="--")



plt.text(510,620,'Main domain', fontsize=10, color="blue")
plt.text(500,580,r'$300km \times 300km  $', fontsize=10, color="blue")
plt.text(550,540,r'$\Delta x, \Delta y$', fontsize=10, color="blue")
plt.text(30,-25,r'$16\Delta x$', fontsize=10, color="blue")
plt.text(200,-25,r'$8\Delta x$', fontsize=10, color="blue")
plt.text(300,-25,r'$4\Delta x$', fontsize=10, color="blue")
plt.text(400,-25,r'$2\Delta x$', fontsize=10, color="blue")
plt.text(570,-25,r'$\Delta x$', fontsize=10, color="blue")

plt.text(1200,90,r'$16\Delta y$', fontsize=10, color="blue")
plt.text(1200,220,r'$8\Delta y$', fontsize=10, color="blue")
plt.text(1200,330,r'$4\Delta y$', fontsize=10, color="blue")
plt.text(1200,410,r'$2\Delta y$', fontsize=10, color="blue")
plt.text(1200,600,r'$\Delta y$', fontsize=10, color="blue")

plt.xlabel("x (km)")
plt.ylabel("y (km)")
plt.title('Gradual horizontal grid')
plt.savefig("./figures/grid.png")

