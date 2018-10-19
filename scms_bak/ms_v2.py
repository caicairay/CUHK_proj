#!/usr/bin/python

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import h5py
import random

def hessian(x):
    """
    Calculate the hessian matrix with finite differences
    Parameters:
       - x : ndarray
    Returns:
       an array of shape (x.dim, x.ndim) + x.shape
       where the array[i, j, ...] corresponds to the second derivative x_ij
    """
    x_grad = np.gradient(x) 
    hessian = np.empty((x.ndim, x.ndim) + x.shape, dtype=x.dtype) 
    for k, grad_k in enumerate(x_grad):
        # iterate over dimensions
        # apply gradient again to every component of the first derivative.
        tmp_grad = np.gradient(grad_k) 
        for l, grad_kl in enumerate(tmp_grad):
            hessian[k, l, :, :] = grad_kl
    return hessian

def neighbour(x,y,lim):
    neighb=[]
    xx=np.arange(x[0]-y,x[0]+y+1)
    yy=np.arange(x[1]-y,x[1]+y+1)
    for i in xx:
        for j in yy:
            dist=np.hypot(x[0]-i,x[1]-j)
            if i>0 and j>0 and i<lim and j<lim and dist<y:
                neighb.append([i,j])
    return neighb

def local_max(pos,rhos,ran):
    n=rhos.shape
    neighb=neighbour(pos,ran,n[0])
    mm=len(neighb)
    rhos_neig=np.repeat(1,mm)
    for (m,ii) in enumerate(neighb):    
        rhos_neig[m]=rhos[ii[0],ii[1]]
    npos=neighb[np.argmax(rhos_neig)]
    shift=npos-pos
    return shift.reshape((2,1))

def plot(t,rhos,mesh,arg):
    plt.figure(t)
    plt.contourf(np.log(rhos),30)
    plt.axis('square')
    for y,x in mesh:
        plt.plot(x,y, 'x', **arg)
    plt.savefig('points'+'{:0>2}'.format(t)+'.jpg',dpi=300,format='jpg',bbox_inches='tight')
    plt.close()
    
def filter(rhos):
    arg={'color':'b', 'markersize':.4} 
    rhos=rhos/np.min(rhos)
    n=rhos.shape
    threshold=np.mean(rhos)*2
    x=np.arange(n[0])
    y=np.arange(n[1])
    X,Y=np.meshgrid(x,y)
    g=np.array(np.gradient(rhos))
    hess=hessian(rhos)
    eps=.2
    ran=5
    filt=np.zeros_like(rhos,dtype=bool) 
    for (i,j),rho in np.ndenumerate(rhos):
        if rhos[i,j]<threshold:
            continue
        shift=local_max(np.array([i,j]),rhos,ran)
        mhess=hess[:,:,i,j]
        w,v=LA.eigh(mhess)
        w2 = w[0]
        if w2>0:
            continue
        vv=v[:,0]
        vv.shape=(2,1)
        pjstep=vv.T.dot(shift)
        if np.abs(pjstep)<=eps:
            filt[i,j]=True
    plt.figure()
    plt.contourf(X,Y,rhos,30)
    for (j,i),x in np.ndenumerate(filt):
        if x == True:
            plt.plot(i,j,'x',color='b',markersize=.6) 
    plt.axis('square')
    plt.savefig('points.jpg',dpi=300,format='jpg',bbox_inches='tight')
    plt.close()
    return filt

def main():
    fnum=17
    path='/data/zitan/240/240-4/'
    #path='/data/zitan/480/480-4/'
    #path='/data/zitan/960/run960-6/'
    flnm=path+'hdfaa.'+'{:03d}'.format(fnum)
    file = h5py.File(flnm,'r')
    rho=file['gas_density']
    rhos=rho[:,:,120]
    point = filter(rhos)

if __name__ == "__main__":
    main()
