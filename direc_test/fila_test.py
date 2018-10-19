import numpy as np
from numpy import linalg as LA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
import toolfunc as tf
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
#from mayavi import mlab
def angle(v1, v2, acute):
# v1 is your firsr vector
# v2 is your second vector
    if (acute == True):
        costheta=np.abs(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    else:

        costheta=np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(costheta)
    return np.rad2deg(angle)

def main():
##____________Control Panel____________##
    method='quiver'
##____________Constants____________##
    colorg=(0,1,0)
    colorr=(1,0,0)
    colorb=(0,0,1)
    colork=(0,0,0)
    x, y, z = np.mgrid[0:4:40j, -4:4:40j, 0:4:40j]

##____________read file_____________##
    file=h5py.File('/data/zitan/960/run960-6/hdfaa.017','r')
    rho=file['gas_density'][:]
    if method=='wholespace':
        rho=rho[280:450,660:900,760:960]
        mlab.figure(bgcolor=(0., 0., 0.), fgcolor=(1., 1., 1.))
        mlab.contour3d(rho,
                       transparent=True,
                       contours=[500,2000,4000,6000,8000,10000],
                       opacity=.5,
                       vmin=500,
                       vmax=10000
                       )
        mlab.outline()
        mlab.axes()
        mlab.title(method)
        mlab.colorbar(orientation='vertical')
    if method=='quiver':
        rho=rho[280:450,660:900,760:960]
        hess=tf.hessian(rho)
        threshold=500
        mesh=rho>threshold
        mesh=np.array([ii for (ii,boo) in np.ndenumerate(mesh) if boo]).astype(int)
        mesh=mesh[::20]
        vv=[]
        for pos in mesh:
            w,v=LA.eigh(hess[(Ellipsis,)+tuple(pos)])
            vv.append(v[:,2])
        vv=np.array(vv)
        ang=[]
        for vec in vv:
            ang.append(angle(vec,[0,0,1],True))
        print 'calculation complete'
        #np.histogram
        plt.hist(ang,bins='auto')
        plt.savefig('hist_orientation.jpg',format='jpg',dpi=250)
        #plt.show()
        #mlab.figure(bgcolor=(0., 0., 0.), fgcolor=(1., 1., 1.))
        #mlab.quiver3d(mesh[:,0],mesh[:,1],mesh[:,2],vv[:,0],vv[:,1],vv[:,2],color=colorr,scale_factor=1.5)
        #mlab.contour3d(rho,
        #               transparent=True,
        #               contours=[500,2000,4000,6000,8000,10000],
        #               opacity=.5,
        #               vmin=500,
        #               vmax=10000
        #               )
        #mlab.outline()
        #mlab.axes()
        #mlab.title(method)
        #mlab.colorbar(orientation='vertical')
    if method=='flow':
        rho=rho[300:340,765:805,830:870]
        hess=tf.hessian(rho)
        threshold=np.mean(rho)
        mask=rho>threshold
        vv=np.zeros((3,40,40,40),dtype=float)
        for (pos,den) in np.ndenumerate(rho):
            w,v=LA.eigh(hess[(Ellipsis,)+tuple(pos)])
            vv[(Ellipsis,)+tuple(pos)]=v[:,2]
        mlab.figure(bgcolor=(0., 0., 0.), fgcolor=(1., 1., 1.))
        mlab.flow(x,y,z,vv[0],vv[1],vv[2])
        mlab.contour3d(x,y,z,rho,transparent=True)
        mlab.outline()
        mlab.axes()
        mlab.title(method)
        mlab.colorbar(orientation='vertical')
    #mlab.show()
if __name__ == "__main__":
    main()
