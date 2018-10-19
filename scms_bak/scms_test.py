#!/usr/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
from numpy import linalg as LA
import matplotlib.pyplot as plt
import h5py
from itertools import product

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

def gaussian(pos,posn,h):
    gaussian=lambda d: 1/np.sqrt(2*np.pi)/h*np.exp(-1/2*(d/h)**2)
    d=np.sqrt(sum((pos-posn)**2))
    return gaussian(d)

def mean_shift(pos,rhos,wholespace,h):
    weight=[gaussian(pos,np.array(ii),h)*den**10 for (ii,den) in np.ndenumerate(rhos)]
    shift=np.average(np.array(wholespace),weights=weight, axis=0)-pos
    return shift[:,np.newaxis]

def plot(t,mesh,rhos):
    fig = plt.figure(t)
    ax = fig.add_subplot(111)#, projection='3d')
    ax.contourf(rhos,30)
    ax.scatter(mesh[:,1], mesh[:,0],s=.1)
    plt.axis('square')
    plt.savefig('points'+'{:0>2}'.format(t)+'.jpg',dpi=300,format='jpg',bbox_inches='tight')
    plt.close()

def scms(rhos):
    wholespace_prod=product(np.arange(rhos.shape[0]),repeat=rhos.ndim)
    wholespace=np.array([ii for ii in wholespace_prod])
    threshold=np.mean(rhos)*3
    rhos_mod=rhos/np.min(rhos)
    mesh=np.zeros_like(rhos,dtype=bool)
    mesh[rhos>threshold]=1
    mesh=np.array([ii for (ii,boo) in np.ndenumerate(mesh) if boo])
    mesh=mesh[::2]
    g=np.array(np.gradient(rhos))
    hess=hessian(rhos)
    mm=mesh.shape[0]
    t=0
    eps=.5
    maxt=10
    h=2.5
    error=np.repeat(1,mm)
    w0=np.repeat(0,mm)
    plot(t,mesh,rhos)
    while np.max(error)>=eps and t<maxt:
        t+=1
        print 'iteration:',t
        for (m,pos) in enumerate(mesh[:,]):
            if error[m]<=eps:
                continue
            shift=mean_shift(pos,rhos_mod,wholespace,h)
            step=np.sqrt(sum(shift**2))
            grad=g[(Ellipsis,)+tuple(pos)]
            mhess=hess[(Ellipsis,)+tuple(pos)]
            w,v=LA.eigh(mhess)
            w0[m]=w[0]>0
            vv=v[:,:-1]
            vvt=vv.dot(vv.T)
            pjshift=vvt.dot(shift)
            pjstep=vv.T.dot(shift)
            error[m]=np.sqrt(np.sum(pjstep**2))
            npos=np.array([x+int(round(y)) for x,y in zip(pos,pjshift)])
            npos=npos.clip(min=0,max=rhos.shape[0]-1)
            mesh[m,]=npos
        plot(t,mesh,rhos)
        print np.mean(error)
    print 'Iteration complete, double checking.'
    for (m,pos) in enumerate(mesh):
        if w0[m] or rhos[tuple(pos)]<threshold:
            mesh[m,]=[0,0]
    mesh=np.unique(mesh,axis=0) 
    print mesh.shape
    print 'Double-check complete, ploting'
    plot(t,mesh,rhos)
    return np.array(mesh)

def main():
    fnum=17
    path='/data/zitan/240/240-4/'
    #path='/data/zitan/480/480-4/'
    #path='/data/zitan/960/run960-6/'
    flnm=path+'hdfaa.'+'{:03d}'.format(fnum)
    file = h5py.File(flnm,'r')
    rho=file['gas_density']
    rhos=rho[:,:,120]
    point = scms(rhos)

if __name__ == "__main__":
    main()

