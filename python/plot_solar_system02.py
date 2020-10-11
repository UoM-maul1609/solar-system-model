import matplotlib
matplotlib.use('agg')
from matplotlib import rc

rc('font',family='serif')
rc('text',usetex = True)

from netCDF4 import Dataset as NetCDFFile

import numpy as np

nc=NetCDFFile('/tmp/output.nc')

import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

objs=('Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto')

# Plot the orbits
plt.ion()

fig=plt.figure()
for i in np.mgrid[1:10:1]:
   ax=plt.subplot(3,3,i)
   ax.plot(nc.variables['time'][:]/365.25/86400.,
       np.sqrt(np.sum(nc.variables['pos'][:,i,0:2]**2,axis=1)),'k.',ms=1)
   ax.set_xlabel('Earth years')
   ax.set_ylabel('%s - Sun distance' % objs[i-1])

nc.close()

mng=plt.get_current_fig_manager()
mng.full_screen_toggle()

fig.savefig('milankovitch.png')
mng.full_screen_toggle()
