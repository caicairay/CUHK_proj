#IMPORTS
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys

### FUNCTIONS ###

### read HDF5 file ###
def read_data(data_num):
    global data_den 
    global data_bdry
    global data_t
    global data_momx
    global data_momy
    global data_gamma
    data_num = str('{:0>4}'.format(data_num))
    file = h5py.File('./g0050_'+data_num+'.h5','r')
    data_den = file['den']
    data_bdry = [file['leftBdry'][0],file['rightBdry'][0]]
    data_t = file['t'][0]
    data_momx = file['momx']
    data_momy = file['momy']
    data_gamma = file['adiGamma']
def plt_den():
    ### set range ###
    x = y = np.linspace(data_bdry[0],data_bdry[1],516)
    X, Y = np.meshgrid(x, y)
    ### draw contour ###
    CS = plt.contourf(X, Y, data_den, 50,cmap=cm.jet)
    plt.xlabel("Time: "+str(data_t)+"Gyr")
    plt.ylabel("kAU")
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


num = int(sys.argv[1])

for i in range(num):
    plt.clf()
    read_data(i)
    plt_den()
    plt_mom()
    plt.pause(0.05)
