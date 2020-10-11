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
#L=5000.

NFFT=np.int(2.**np.ceil(np.log(L)/np.log(2.)))
dt=(nc.variables['time'][1]-nc.variables['time'][0]) / (86400.0*365.25)
Fs=1./dt
f = Fs/2.*np.linspace(0.,1.,np.int(NFFT/2+1))

masses=np.array([1.327124400e11,22032.09, 324858.63, \
                 398600.440, 42828.3, 126686511, \
                 37931207.8 , 5793966, 6835107, 872.4])
totalmass=np.sum(masses)

n_bodies=10
plt.ion()
fig=plt.figure()
ax=plt.axes()
sum1=0.
x=1./totalmass*masses[0]*nc.variables['pos'][:,0,0]
y=1./totalmass*masses[0]*nc.variables['pos'][:,0,1]
z=1./totalmass*masses[0]*nc.variables['pos'][:,0,2]

for i in np.mgrid[1:n_bodies:1]:
    x=x+1/totalmass*masses[i]*nc.variables['pos'][:,i,0]
    y=y+1/totalmass*masses[i]*nc.variables['pos'][:,i,1]
    z=z+1/totalmass*masses[i]*nc.variables['pos'][:,i,2]

y1=np.sqrt(x**2+y**2+z**2)


y=np.fft.fft(y1,int(NFFT))/L
ax.plot(1./f,(2.*np.abs(y[1:int(NFFT/2+2)]) ))

ax.semilogx()
ax.semilogy()
plt.xlabel('time period (years)')
plt.ylabel('power')
plt.title('Sun''s Wobble')


nc.close()

fig.savefig('fourier_sun.png')
