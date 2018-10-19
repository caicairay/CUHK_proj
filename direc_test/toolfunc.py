import numpy as np
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

def mean_shift(pos,rhos,ran):
    n=rhos.shape
    neighb=neighbour(pos,ran,n[0])
    weight=[rhos[ii[0],ii[1]]**3 for ii in neighb]
    shiftx=np.average([ii[0] for ii in neighb],weights=weight)-pos[0]
    shifty=np.average([ii[1] for ii in neighb],weights=weight)-pos[1]
    return [shiftx,shifty]

def local_max(pos,rhos,ran):
    n=rhos.shape
    neighb=neighbour(pos,ran,n[0])
    mm=len(neighb)
    rhos_neig=np.repeat(1.,mm)
    for (m,ii) in enumerate(neighb):    
        rhos_neig[m]=rhos[tuple(ii)]
    npos=neighb[np.argmax(rhos_neig)]
    shift=np.array(npos)-np.array(pos)
    return shift
