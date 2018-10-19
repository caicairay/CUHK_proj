import numpy as np
import scipy.integrate as integrate
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import h5py

AU = 1.496e13 # [cm]
Msol = 1.99e33 # [g]
year = 3.156e7 # [sec]

def read_data(data_num):
    global ndom,nmesh,x,y,xx,yy,r,velx,vely,vr,vt,data_den,sgfx,sgfy,sgfr,sgft
    data_num = str('{:0>4}'.format(data_num))
    file = h5py.File('./data/g0050_'+data_num+'.h5','r')
    ndom = int( file['nMesh'][0])
    nmesh = ndom +4
    data_den= file['den'][2:ndom+2,2:ndom+2]
    data_bdry=[file['leftBdry'][:],file['rightBdry'][:]]
    data_t = file['t'][:]
    data_momx = file['momx'][2:ndom+2,2:ndom+2]
    data_momy = file['momy'][2:ndom+2,2:ndom+2]
    sgfx = file['sgfx'][:]
    sgfy = file['sgfy'][:]
    data_pos = [file['xc1'][2:ndom+2],file['xc2'][2:ndom+2]]
    x,y = data_pos
    xx,yy = np.meshgrid(x,y)
    r = np.sqrt(xx**2+yy**2)
    velx=data_momx/data_den
    vely=data_momy/data_den
    theta = np.arctan2(yy,xx)
    vr = np.cos(theta)*velx+np.sin(theta)*vely
    vt = -np.sin(theta)*velx+np.cos(theta)*vely
    sgfr = np.cos(theta)*sgfx+np.sin(theta)*sgfy
    sgft = -np.sin(theta)*sgfx+np.cos(theta)*sgfy

def plt_vr():
    plt.figure()
    plt.plot(r[ndom/2][:ndom/2]/AU,vr[ndom/2][:ndom/2],'.')
    #plt.plot(r[ndom/2][ndom/2:],vr[ndom/2][ndom/2:],'.')
    plt.xlabel("r")
    plt.ylabel(r"$v_r$")
    plt.savefig('../readdata/vr.png')

def plt_vt():
    plt.figure()
    plt.plot(r[ndom/2][:ndom/2]/AU,vt[ndom/2][:ndom/2],'.')
    #plt.plot(r[ndom/2][ndom/2:],vt[ndom/2][ndom/2:],'.')
    plt.xlabel("r")
    plt.ylabel(r"$v_\theta$")
    plt.savefig('../readdata/vt.png')

def plt_den():
    plt.figure()
    plt.plot(r[ndom/2][:ndom/2]/AU,np.log10(data_den[ndom/2][:ndom/2]),'.')
    #plt.plot(r[ndom/2][ndom/2:],data_den[ndom/2][ndom/2:],'.')
    plt.xlabel("r")
    plt.ylabel(r"$\rho$")
    plt.savefig('../readdata/rho.png')

def plt_grav():
    plt.figure()
    plt.plot(r[ndom/2][:ndom/2]/AU,sgfr[ndom/2][:ndom/2],'.')
    plt.xlabel("r")
    plt.ylabel("gravr")
    plt.figure()
    plt.plot(r[ndom/2][:ndom/2],sgft[ndom/2][:ndom/2],'.')
    plt.xlabel("r")
    plt.ylabel("gravt")
    plt.savefig('../readdata/grav.png')

if __name__ == "__main__":
    data_num = sys.argv[1]
    read_data(data_num)
    plt_vr()
    plt_vt()
    plt_den()
    plt_grav()
