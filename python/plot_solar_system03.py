import matplotlib
matplotlib.use('agg')
from matplotlib import rc

rc('font',family='serif')
rc('text',usetex = True)

from netCDF4 import Dataset as NetCDFFile

import numpy as np

nc=NetCDFFile('/tmp/output.nc')

import matplotlib.pyplot as plt




# Do Fourier analysis:
L=50000.
#L=600000.
#L=40000.

NFFT=np.int(2.**np.ceil(np.log(L)/np.log(2.)))
dt=(nc.variables['time'][1]-nc.variables['time'][0]) / (86400.0*365.25)
Fs=1./dt
f = Fs/2.*np.linspace(0.,1.,np.int(NFFT/2+1))

n_bodies=10
plt.ion()
fig=plt.figure()
ax=plt.axes()
for i in np.mgrid[1:n_bodies:1]:
    y1=np.sqrt(np.sum(nc.variables['pos'][:,i,0:2].astype(float)**2,axis=1))
    y=np.fft.fft(y1,int(NFFT))/L
    ax.plot(1./f,(2.*np.abs(y[1:int(NFFT/2+2)]) ))
#   ax.plot(f,(2.*np.abs(y[1:NFFT/2+2]) ))

ax.semilogx()
ax.semilogy()
ax.legend(('Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto'))
plt.xlim((10,30))
plt.xlabel('time period (years)')
plt.ylabel('power')
"""
fig=plt.figure()
ax=plt.axes()
for i in np.mgrid[1:n_bodies:1]:
   y=np.fft.fft(np.sqrt((nc.variables['pos'][:,i,0]**2)),int(NFFT))/L
   ax.plot(1./f,(2.*np.abs(y[1:int(NFFT/2+2)]) ))

ax.semilogx()
ax.semilogy()
ax.legend(('Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto'))



"""

nc.close()

fig.savefig('fourier.png')
