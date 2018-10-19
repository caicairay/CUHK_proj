#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import h5py
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

def hessian(x):
    """
    Calculate the hessian matrix with finite differences
    Parameters:
       - x : ndarray
    Returns:
       an array of shape (x.ndim, x.ndim) + x.shape
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

def angle(v1, v2, acute):
# v1 is your firsr vector
# v2 is your second vector
    if (acute == True):
        costheta=np.abs(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    else:

        costheta=np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(costheta)
    return np.rad2deg(angle)

def filter(rhos):
    n=rhos.shape
    threshold=np.mean(rhos)*1.5
    g=np.array(np.gradient(rhos))
    hess=hessian(rhos)
    eps=15
    filt=np.zeros_like(rhos,dtype=bool) 
    for (pos,rho) in np.ndenumerate(rhos): 
        grad=g[(Ellipsis,)+pos]
        mhess=hess[(Ellipsis,)+pos]
        w,v=LA.eigh(mhess)
        w0 = w[0]
        vv=v[:,-1]
        deg=angle(grad,vv,False)
        deg=np.nan_to_num(deg)
        if deg<=eps and w0<0 and rhos[pos]>threshold:
            filt[pos]=True
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #plt.contourf(rhos,30)
    xs=[]
    ys=[]
    zs=[]
    for (j,i,k),x in np.ndenumerate(filt):
        if x == True:
            xs.append(i)
            ys.append(j)
            zs.append(k)
    ax.scatter(xs, ys, zs,s=.5)
    plt.savefig('filt.jpg',dpi=300,format='jpg',bbox_inches='tight')
    plt.close()
    return filt

def main():
    fnum=17
    #path='/data/zitan/240/240-4/'
    #path='/data/zitan/480/480-4/'
    path='/data/zitan/960/run960-6/'
    flnm=path+'hdfaa.'+'{:03d}'.format(fnum)
    file = h5py.File(flnm,'r')
    rho=file['gas_density']
    x1=300;
    x2=339;
    y1=765;
    y2=804;
    z1=830;
    z2=869;
    rhos=rho[x1:x2,y1:y2,z1:z2]
    point = filter(rhos)

if __name__ == "__main__":
    main()
