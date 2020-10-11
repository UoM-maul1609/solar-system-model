from matplotlib import rc

rc('font',family='serif')
rc('text',usetex = True)

from netCDF4 import Dataset as NetCDFFile

import numpy as np

nc=NetCDFFile('/tmp/output.nc')

import matplotlib.pyplot as plt

G=6.67e-11
m=1.989e30

Gm=np.array([G*m/1e9, 
   22032.09, 324858.63, 398600.440,
    42828.3, 126686511, 37931207.8 ,
    5793966, 6835107, 872.4])*1.e9

TIME=nc.variables['time'][:]
(r,)=np.shape(TIME)

x=nc.variables['pos'][:,:,0] #-np.reshape(np.repeat(nc.variables['POS'][:,0,0],10),[r,10])
y=nc.variables['pos'][:,:,1] #-np.reshape(np.repeat(nc.variables['POS'][:,0,1],10),[r,10])
z=nc.variables['pos'][:,:,2] #-np.reshape(np.repeat(nc.variables['POS'][:,0,2],10),[r,10])

ux=nc.variables['vel'][:,:,0] #-np.reshape(np.repeat(nc.variables['VEL'][:,0,0],10),[r,10])
uy=nc.variables['vel'][:,:,1] #-np.reshape(np.repeat(nc.variables['VEL'][:,0,1],10),[r,10])
uz=nc.variables['vel'][:,:,2] #-np.reshape(np.repeat(nc.variables['VEL'][:,0,2],10),[r,10])


lx=(y*uz-z*uy)*Gm/G
ly=(z*ux-x*uz)*Gm/G
lz=(x*uy-y*ux)*Gm/G
lmag=np.sqrt(lx*lx+ly*ly+lz*lz)
nc.close()

