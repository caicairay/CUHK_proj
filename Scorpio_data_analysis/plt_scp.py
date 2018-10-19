#imports
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
### read HDF5 file ###
data_num = sys.argv[1]
data_num = str('{:0>4}'.format(data_num))
file = h5py.File('./g0050_'+data_num+'.h5','r')
data_den= file['den']
data_bdry= [file['leftBdry'][0],file['rightBdry'][0]]
data_t = file['t'][0]
data_momx = file['momx']
data_momy = file['momy']
print data_den[1][1]
def plt_den():
    #set range
    x = y = np.linspace(data_bdry[0],data_bdry[1],516)
    X, Y = np.meshgrid(x, y)
    #draw contour
    CS = plt.contourf(X, Y, data_den, 50,cmap=cm.jet)
    plt.xlabel("Time: "+str(data_t)+"Gyr")
    plt.ylabel("kpc")
    plt.colorbar()

def plt_mom():
    x = y = np.linspace(data_bdry[0],data_bdry[1],516)
    X, Y = np.meshgrid(x, y)
    M = np.hypot(data_momx, data_momy)
    Q = plt.quiver(X[::10,::10], Y[::10,::10], data_momx[::10,::10],
                   data_momy[::10,::10], M[::10,::10], pivot='tail',
                   units='width')
    plt.colorbar()
    qk = plt.quiverkey(Q, 0.9, 0.9, 2, r'$\frac{Msun}{AU^2}\times\frac{km}{s}$', labelpos='E',
                   coordinates='figure')

print file.items()
plt_den()
plt_mom()
plt.show()
file.close()
