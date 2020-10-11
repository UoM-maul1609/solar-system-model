from matplotlib import rc

rc('font',family='serif')
rc('text',usetex = True)

from netCDF4 import Dataset as NetCDFFile

import numpy as np

nc1=NetCDFFile('output03.nc')
nc2=NetCDFFile('output04.nc')
nc3=NetCDFFile('output05.nc')

import matplotlib.pyplot as plt

from scipy.signal import find_peaks_cwt
ii=1

plt.ion()
fig=plt.figure()
ax=plt.axes()

dist1=np.sqrt(np.sum((nc1.variables['POS'][:,ii,0:2]-nc1.variables['POS'][:,0,0:2])**2,axis=1))
dist2=np.sqrt(np.sum((nc2.variables['POS'][:,ii,0:2]-nc2.variables['POS'][:,0,0:2])**2,axis=1))
dist3=np.sqrt(np.sum((nc3.variables['POS'][:,ii,0:2]-nc3.variables['POS'][:,0,0:2])**2,axis=1))

TIME=nc1.variables['TIME'][:]/365.25/86400

#
d1=np.diff(dist1,n=1)
i1=d1[1:len(d1)]*d1[0:len(d1)-1]
ind1=(i1<0.)
TIME1=TIME[ind1];dist11=dist1[ind1]
indexes1=np.arange(1,len(TIME1),2)
#
d2=np.diff(dist2,n=1)
i2=d2[1:len(d2)]*d2[0:len(d2)-1]
ind2=(i2<0.)
TIME2=TIME[ind2];dist22=dist2[ind2]
indexes2=np.arange(1,len(TIME2),2)
#
d3=np.diff(dist3,n=1)
i3=d3[1:len(d3)]*d3[0:len(d3)-1]
ind3=(i3<0.)
TIME3=TIME[ind3];dist33=dist3[ind3]
indexes3=np.arange(1,len(TIME3),2)

#indexes1=find_peaks_cwt(dist1[1:50000],np.arange(0.1,1.))
#indexes2=find_peaks_cwt(dist2[1:50000],np.arange(0.1,1.))
#indexes3=find_peaks_cwt(dist3[1:50000],np.arange(0.1,1.))

y1=nc1.variables['POS'][:,ii,1]-nc1.variables['POS'][:,0,1];y1=y1[ind1]
x1=nc1.variables['POS'][:,ii,0]-nc1.variables['POS'][:,0,0];x1=x1[ind1]

y2=nc2.variables['POS'][:,ii,1]-nc2.variables['POS'][:,0,1];y2=y2[ind2]
x2=nc2.variables['POS'][:,ii,0]-nc2.variables['POS'][:,0,0];x2=x2[ind2]

y3=nc3.variables['POS'][:,ii,1]-nc3.variables['POS'][:,0,1];y3a=y3;y3=y3[ind3]
x3=nc3.variables['POS'][:,ii,0]-nc3.variables['POS'][:,0,0];x3a=x3;x3=x3[ind3]

plt.plot(TIME1[indexes1],180./np.pi*np.arctan(y1[indexes1]/x1[indexes1]),'r')

plt.plot(TIME2[indexes2],180./np.pi*np.arctan(y2[indexes2]/x2[indexes2]),'g')

plt.plot(TIME3[indexes3],180./np.pi*np.arctan(y3[indexes3]/x3[indexes3]),'b')

p1=np.polyfit(TIME1[indexes1[:]],180./np.pi*np.arctan(y1[indexes1[:]]/x1[indexes1[:]]),1)
p2=np.polyfit(TIME2[indexes2[:]],180./np.pi*np.arctan(y2[indexes2[:]]/x2[indexes2[:]]),1)
p3=np.polyfit(TIME3[indexes3[:]],180./np.pi*np.arctan(y3[indexes3[:]]/x3[indexes3[:]]),1)


print p1[0]*3600.,p2[0]*3600.,p3[0]*3600.
nc1.close()
nc2.close()
nc3.close()

