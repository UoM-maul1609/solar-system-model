import os
import getpass


import matplotlib
matplotlib.use('agg')
from matplotlib import rc

rc('font',family='serif')
rc('text',usetex = True)

from netCDF4 import Dataset as NetCDFFile

import numpy as np


username=getpass.getuser()

nc=NetCDFFile('/tmp/' + username + '/output.nc')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


r,=np.shape(nc.variables['time'][:])

# Plot the orbits, every 10th point
stride=1
plt.ion()
fig=plt.figure()
ax=Axes3D(fig,auto_add_to_figure=False)
fig.add_axes(ax)
for i in np.mgrid[1:10:1]:
   ax.plot(nc.variables['pos'][1:r:stride,i,0],nc.variables['pos'][1:r:stride,i,1],
         nc.variables['pos'][1:r:stride,i,2],'k.',ms=0.1)


ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title('Orbits of the planets')
nc.close()

mng=plt.get_current_fig_manager()
mng.full_screen_toggle()

fig.savefig('/tmp/' + username + '/orbits.png', dpi=300)

mng.full_screen_toggle()

