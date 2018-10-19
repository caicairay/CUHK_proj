import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import h5py
import toolfunc as tf
def main():
    rhos=np.loadtxt('../data_hub/ring.txt')
    x=np.linspace(0,1,rhos.shape[0])
    y=np.linspace(0,1,rhos.shape[1])
    X,Y=np.meshgrid(x,y)

    hess=tf.hessian(rhos)

    vv0=np.zeros((2,rhos.shape[0],rhos.shape[1]))
    vv1=np.zeros((2,rhos.shape[0],rhos.shape[1]))
    pjshift0=np.zeros((2,rhos.shape[0],rhos.shape[1]))
    pjshift1=np.zeros((2,rhos.shape[0],rhos.shape[1]))
    g=np.zeros((2,rhos.shape[0],rhos.shape[1]))

    mask=np.zeros((rhos.shape[0],rhos.shape[1]),dtype=bool)

    for (idx,den) in np.ndenumerate(rhos):
        grad=tf.mean_shift(idx,rhos,5)
        g[(Ellipsis,)+idx]=grad
        w,v=LA.eigh(hess[(Ellipsis,)+idx])
        mask[idx]=w[0]<0 
        vv0[(Ellipsis,)+idx]=v[:,0]
        vv1[(Ellipsis,)+idx]=v[:,1]
        vv=v[:,0].reshape((2,1))
        vvt=vv.dot(vv.T)
        pjshift0[(Ellipsis,)+idx]=vvt.dot(grad)
        vv=v[:,1].reshape((2,1))
        vvt=vv.dot(vv.T)
        pjshift1[(Ellipsis,)+idx]=vvt.dot(grad)
    pjstep0=np.hypot(pjshift0[0],pjshift0[1])
    pjstep1=np.hypot(pjshift1[0],pjshift1[1])

    U_v0=vv0[0]
    V_v0=vv0[1]
    U_v1=vv1[0]
    V_v1=vv1[1]
    U_g=g[0]
    V_g=g[1]
    #U_pj0=pjshift0[0]/pjstep0
    #V_pj0=pjshift0[1]/pjstep0
    U_pj0=pjshift0[0]
    V_pj0=pjshift0[1]
    #U_pj1=pjshift1[0]/pjstep1
    #V_pj1=pjshift1[1]/pjstep1
    U_pj1=pjshift1[0]
    V_pj1=pjshift1[1]

    #plt.contourf(X,Y,rhos,50)
    plt.contourf(X[::10,::10],Y[::10,::10],rhos[::10,::10],50)

    plt.quiver(X[::10,::10],Y[::10,::10],U_v0[::10,::10],V_v0[::10,::10],color='r')
    plt.quiver(X[::10,::10],Y[::10,::10],U_v1[::10,::10],V_v1[::10,::10],color='b')
    plt.quiver(X[::10,::10],Y[::10,::10],U_g[::10,::10],V_g[::10,::10])

    #plt.quiver(X[::10,::10][mask[::10,::10]],Y[::10,::10][mask[::10,::10]],U_v0[::10,::10][mask[::10,::10]],V_v0[::10,::10][mask[::10,::10]],color='r')
    #plt.quiver(X[::10,::10][mask[::10,::10]],Y[::10,::10][mask[::10,::10]],U_v1[::10,::10][mask[::10,::10]],V_v1[::10,::10][mask[::10,::10]],color='b')
    #plt.quiver(X[::10,::10][mask[::10,::10]],Y[::10,::10][mask[::10,::10]],U_g[::10,::10][mask[::10,::10]],V_g[::10,::10][mask[::10,::10]])
    #plt.quiver(X[::10,::10],Y[::10,::10],U_g[::10,::10],V_g[::10,::10],angles='xy')

    #plt.quiver(X[::10,::10][mask[::10,::10]],Y[::10,::10][mask[::10,::10]],U_pj0[::10,::10][mask[::10,::10]],V_pj0[::10,::10][mask[::10,::10]],color='cyan')
    #plt.quiver(X[::10,::10][mask[::10,::10]],Y[::10,::10][mask[::10,::10]],U_pj1[::10,::10][mask[::10,::10]],V_pj1[::10,::10][mask[::10,::10]],color='g')

    plt.grid(b=True)
    plt.axis('square')
    plt.show()
    #plt.savefig('dirc_test.jpg',format='jpg',dpi=250)
if __name__ == "__main__":
    main()
