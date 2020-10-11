from matplotlib import rc

rc('font',family='serif')
rc('text',usetex = True)

from netCDF4 import Dataset as NetCDFFile

import numpy as np

nc=NetCDFFile('output02.nc')

import matplotlib.pyplot as plt


G=6.67e-11
m=np.array([1.989e30, 328.5e21, 4.867e24, 100e30])

Gm=np.array([G*m[1]/1e9, 22032.09, 324858.63, 398600.440,
    42828.3, 126686511, 37931207.8,
    5793966, 6835107, 872.4, G*5e30/1e9])*1e9
Rsun=6.955e8
dR=Rsun

r,c,p=np.shape(nc.variables['POS'][:][:])

tide_x=np.zeros((r,c-1))
tide_y=np.zeros_like(tide_x)
tide_z=np.zeros_like(tide_x)

for i in np.mgrid[1:10:1]:
#   if (i==5 or i==6 ):
   if (i >=1):
      # tidal forces in the x,y,z directions
      tide_x[:,i-1]=dR*2.*Gm[i] / (np.sum(nc.variables['POS'][:,i,0:2]**2,axis=1))**2 * nc.variables['POS'][:,i,0]
      tide_y[:,i-1]=dR*2.*Gm[i] / (np.sum(nc.variables['POS'][:,i,0:2]**2,axis=1))**2 * nc.variables['POS'][:,i,1]
      tide_z[:,i-1]=dR*2.*Gm[i] / (np.sum(nc.variables['POS'][:,i,0:2]**2,axis=1))**2 * nc.variables['POS'][:,i,2]

tide_mag=np.sqrt(np.sum(tide_x,axis=1)**2+np.sum(tide_y,axis=1)**2+
                 np.sum(tide_z,axis=1)**2)


# Plot it all out:
plt.ion()
plt.figure()
TIME=nc.variables['TIME'][1:50000]/86400/365.25
plt.plot(TIME,tide_mag[1:50000],'k',lw=1)

from scipy.signal import filtfilt
import scipy.signal as signal

b, a = signal.butter(8, 0.125)


ydat=filtfilt(b,a,tide_mag[1:50000])
plt.plot(TIME,ydat,'r',lw=1)

plt.xlim((0,100))
plt.xlabel('Years')
plt.ylabel('Magnitude of tidal force on sun')
plt.legend(('Raw','Smoothed'))

plt.savefig('time_series_of_tides.png')



# Do a Fourier analysis of the signal 
L=6000.
L=30000.
L=100000.

NFFT=2.**np.ceil(np.log(L)/np.log(2.))
dt=TIME[1]-TIME[0]
Fs=1./dt
f = Fs/2.*np.linspace(0,1,NFFT/2.+1.)

fig=plt.figure()
ax=plt.axes()
#y=np.fft.fft(ydat,int(NFFT))/L
y=np.fft.fft(tide_mag,int(NFFT))/L
ax.plot(1./f,(2.*np.abs(y[1:NFFT/2+2]) ))

#ax.semilogx()
ax.semilogy()
ax.set_xlim((0.,40.))

plt.plot((11.,11.),plt.ylim(),'k')
plt.plot((20.,20.),plt.ylim(),'r')
plt.plot((22.,22.),plt.ylim(),'k')

plt.xlabel('Years')
plt.ylabel('Power')


fig.savefig('fourier_of_tides.png')


# Allow signals with frequency lower than 1/10 yr-1 pass
# http://stackoverflow.com/questions/19122157/fft-bandpass-filter-in-python
import numpy as np
from scipy.fftpack import rfft, irfft, fftfreq
TIME=nc.variables['TIME'][:]/86400/365.25

W = fftfreq(tide_mag.size, d=TIME[1]-TIME[0])
f_signal = rfft(tide_mag)

# If our original signal time was in seconds, this is now in Hz    
cut_f_signal = f_signal.copy()
#cut_f_signal[(1./W<2.)] = 0
cut_f_signal[(1./W<8.)] = 0

cut_signal = irfft(cut_f_signal)
plt.figure()
#plt.plot(TIME[0:len(cut_signal)+1],cut_signal)
plt.plot(2014.-TIME[0:len(cut_signal)+1],cut_signal)


#plt.plot(2014-nc.variables['TIME'][:]/365.25/86400,np.sum((nc.variables['POS'][:,6,0:3]-nc.variables['POS'][:,5,0:3])**2,axis=1) )
nc.close()

