#!/usr/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
from numpy import linalg as LA
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
from itertools import product
from sklearn import preprocessing

##____________SCMS____________##
def scms(rhos,plot_path,maxt=10,threshold=None,plot_switch=False,win_mode='gaussian',map_mode='expo',win_para=2.5,map_para=0.3):
    """
    Parameters:
        - rhos: density:square array
        - maxt: the maximum iteration number
        - threshold: the threshold of the density
        - win_mode & win_para:
            - gaussian(h)
            - flat(ran)
        - map_mode & map_para:
            - polynomial(gamma)
            - sigmoid(threshold)
    Returns:
        - points: n-D points, represent the positions of the ridges
    """
    ##___________PlotFunction__________##
    def plot(t,mesh,rhos,path):
        fig = plt.figure(t)
        ax = fig.add_subplot(111)#, projection='3d')
        ax.contourf(rhos,30)
        #ax.scatter(mesh[:,1],mesh[:,0],s=.1,edgecolors='none')
        ax.plot(mesh[:,1],mesh[:,0],'x',color='r',markersize=0.2)
        plt.axis('square')
        plt.savefig(path+'/points'+'{:0>2}'.format(t)+'.jpg',dpi=250,format='jpg',bbox_inches='tight')
        plt.close()
    ##____________HessianCalculation_____________##
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
            """
            iterate over dimensions
            apply gradient again to every component of the first derivative.
            """
            tmp_grad = np.gradient(grad_k) 
            for l, grad_kl in enumerate(tmp_grad):
                hessian[k, l, :, :] = grad_kl
        return hessian
    
    ##___________WeightCalculation___________##
    def weight_calculation(pos,rhos,win_mode,map_mode,win_para,map_para):
        """
        Window:Gaussian(h) or Flat(ran) or Both([h,ran])
        Map_func:Polynomial(gamma) or Sigmoid(threshold) or Exponential(sigma)
        """
        ##___________WindowFunctions___________##
        def gaussian(pos,posn,h):
            gaussian=lambda d: 1/np.sqrt(2*np.pi)/h*np.exp(-1/2*(d/h)**2)
            d=np.sqrt(sum((pos-posn)**2))
            return gaussian(d)
        def flat(pos,posn,ran):
            dist=np.sqrt(sum((pos-posn)**2))
            if dist<ran:
                return 1
            else:
                return 0
        def comb_win(pos,posn,hran):
            return gaussian(pos,posn,hran[0])*flat(pos,posn,hran[1])
        ##___________MapFunctions___________##
        def poly(den,gamma):
            return den**gamma
        def sigmoid(den,threshold):
            return 1./(1+np.exp(-(den-threshold)))
        def expo(den,sigma):
            return np.exp(den/sigma**2)
        windict={
            'flat':flat,
            'gaussian':gaussian,
            'both':comb_win
                 }
        mapdict={
            'poly':poly,
            'sigmoid':sigmoid,
            'expo':expo
                }
        weight=[windict[win_mode](pos,posn,win_para)*mapdict[map_mode](den,map_para) for (posn,den) in np.ndenumerate(rhos)]
        return weight
    
    ##___________MeanShiftImplement___________##
    def mean_shift(wholespace,weight):
        shift=np.average(np.array(wholespace),weights=weight, axis=0)-pos
        return shift[:,np.newaxis]

##_____________PreWork______________#
    if threshold==None:
        threshold=np.mean(rhos)
    wholespace_prod=product(np.arange(rhos.shape[0]),repeat=rhos.ndim)#
    wholespace=np.array([ii for ii in wholespace_prod])               #Generate positions for wholespace
    rhos_norm=preprocessing.normalize(rhos,norm='max')                 #Normalizing data for Mean Shift
    mesh=rhos>threshold                                                 #
    mesh=np.array([ii for (ii,boo) in np.ndenumerate(mesh) if boo])     #
    mesh=mesh[::2]                                                      #Generate start points
    g=np.array(np.gradient(rhos))
    hess=hessian(rhos)
##______________SCMSLOOP_______________##
    mm=mesh.shape[0]
    t=0
    eps=.5
    error=np.repeat(1,mm)
    w0=np.repeat(0,mm)
    if plot_switch==True:
        plot(t,mesh,rhos,plot_path)
    while np.max(error)>=eps and t<maxt:
        t+=1
        print 'iteration:',t
        for (m,pos) in enumerate(mesh[:,]):
            if error[m]<=eps:
                continue
            weight=weight_calculation(pos,rhos_norm,win_mode,map_mode,win_para,map_para)
            shift=mean_shift(wholespace,weight)
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
        if plot_switch==True:
            plot(t,mesh,rhos,plot_path)
        print 'Mean error:'+str(np.mean(error))
##_____________DoubleCheck______________##
    print 'Iteration complete, double checking.'
    for (m,pos) in enumerate(mesh):
        if w0[m] or rhos[tuple(pos)]<threshold:
            mesh[m,]=None
    mesh=[x for x in mesh if x is not None]
    mesh=np.unique(mesh,axis=0) 
##______________WriteGreating____________##
    print 'Double-check complete,'+str(mesh.shape[0])+'/'+str(mm)+' points left'
    if plot_switch==True:
        print 'Plotting'
        plot(99,mesh,rhos,plot_path)
    return np.array(mesh)


import h5py
def main():
##________ReadData________##
    fnum=17
    path='/data/zitan/240/240-4/'
    #path='/data/zitan/480/480-4/'
    #path='/data/zitan/960/run960-6/'
    flnm=path+'hdfaa.'+'{:03d}'.format(fnum)
    file = h5py.File(flnm,'r')
    rho=file['gas_density']
##________Slicing_________##
    rhos=rho[:,:,120]
##________ControlPanel_________##
    args={
        'plot_path':'./test',
        'threshold':np.mean(rhos)*3,
        'plot_switch':True,
        'win_mode':'both',
        'win_para':[0.5,5]
        }
##________FilamentDetection_______##
    point = scms(rhos,**args)
    #plot(99,point,rhos)

if __name__ == "__main__":
    main()
