#!/usr/bin/python

import matplotlib
matplotlib.use('Agg')
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import h5py

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

def gaussian(pos,posn,h):
    gaussian=lambda d: 1/np.sqrt(2*np.pi)/h*np.exp(-1/2*(d/h)**2)
    d=np.sqrt(sum((pos-posn)**2))
    return gaussian(d)

def sigmoid(x):
    return 1/(1+np.exp(-x))

def mean_shift(pos,rhos,ran,h):
    n=rhos.shape
    neighb=neighbour(pos,ran,n[0])
    #weight=[(rhos[ii[0],ii[1]]**50)*gaussian(pos,ii,h) for ii in neighb] #
    weight=[(rhos[ii[0],ii[1]]**50) for ii in neighb] #
    shiftx=np.average([ii[0] for ii in neighb],weights=weight)-pos[0]
    shifty=np.average([ii[1] for ii in neighb],weights=weight)-pos[1]
    return np.array([[shiftx],
                     [shifty]])

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

def plot(t,mesh,rhos):
    fig = plt.figure(t)
    ax = fig.add_subplot(111)#, projection='3d')
    ax.contourf(rhos,30)
    ax.scatter(mesh[:,1], mesh[:,0],s=.1)
    plt.axis('square')
    plt.savefig('points'+'{:0>2}'.format(t)+'.jpg',dpi=300,format='jpg',bbox_inches='tight')
    plt.close()

def scms(rhos):
    #arg={'color':'b', 'markersize':.4} 
    #rhos_mod=rhos-np.mean(rhos)
    rhos_mod=rhos/np.min(rhos)
    n=rhos.shape
    threshold=np.mean(rhos)
    mesh=[]
    for i in np.arange(0,n[0]-1,1):
        for j in np.arange(0,n[0]-1,1):
            if rhos[i,j]>threshold:
                mesh.append([i,j])
    mesh=np.array(mesh)
    g=np.array(np.gradient(rhos))
    hess=hessian(rhos)
    mm=mesh.shape[0]
    t=0
    eps=.5
    maxt=10
    ran=np.repeat(5,mm)
    h=np.repeat(15,mm)
    error=np.repeat(1,mm)
    ww=np.repeat(1,mm)
    plot(t,mesh,rhos)
    points=[]
    while max(error)>=eps and t<maxt:
        t+=1
        print 'iteration:',t
        for (m,pos) in enumerate(mesh[:,]):
            if error[m]<=eps:
                continue
            shift=mean_shift(pos,rhos_mod,ran[m],h[m])
            #shift=local_max(pos,rhos,ran[m])
            round_pos=np.array([int(round(ii)) for ii in pos])
            mhess=hess[(Ellipsis,)+tuple(round_pos)]
            w,v=LA.eigh(mhess)
            ww[m]=w[0]
            vv=v[:,0]
            vv.shape=(2,1)
            vvt=vv.dot(vv.T)
            pjshift=vvt.dot(shift)
            pjstep=vv.T.dot(shift)
            error[m]=np.abs(pjstep)
            npos=np.array([x+y for x,y in zip(pos,pjshift)])
            npos=npos.clip(min=0,max=n[0]-1).reshape([2,])
            mesh[m,]=npos
        #plot(t,mesh,rhos)
        print np.mean(error)
    print 'Iteration complete, double checking.'
    for (m,pos) in enumerate(mesh[:,]):
        if ww[m]>0 or rhos[pos[0],pos[1]]<threshold:
            mesh[m,]=[0,0]
    mesh=np.unique(mesh,axis=0) 
    print 'Double-check complete',mesh.shape[0],'/',mm, 'points left'
    print 'Ploting'
    plot(99,mesh,rhos)
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

