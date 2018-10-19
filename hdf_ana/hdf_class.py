import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.axes
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#CONSTANTS
AU = 1.496e13 # [cm]
Msol = 1.99e33 # [g]
year = 3.156e7 # [sec]

#__________#
l_unit=AU
t_unit=year
m_unit=Msol
#^^^^^^^^^^#


def custom_cm():
    cmap_name='my_cm'
    #color_list=[(0,0,0),(0.1,0.2,0.9),(1,1,1),(0.9,0.1,0.2)]
    color_list=[(0,(0,0,0)),(0.2,(0.1,0.2,0.8)),(0.5,(1,1,1)),(1,(0.8,0.1,0.2))]
    nbins=100
    my_cm = LinearSegmentedColormap.from_list(
            cmap_name, color_list, N=nbins)
    return my_cm
def colorbar_right(ax):
    axins = inset_axes(ax,width="3%",  # width = 10% of parent_bbox width
                       height="50%",  # height : 50%
                       loc=4,
                       bbox_to_anchor=(0., 0., 1, 1),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
    return axins
def colorbar_top(ax):
    axins = inset_axes(ax,
                       width="50%",  # width = 10% of parent_bbox width
                       height="2%",  # height : 50%
                       loc=2,
                       bbox_to_anchor=(0., 0.03, 1., 1.),
                       bbox_transform=ax.transAxes,
                       borderpad=0,
                       )
    return axins

class HDF5_DATA():
    def __init__(self,data_num=None):
        data_num = str('{:0>4}'.format(data_num))
        file = h5py.File('../data/g0050_'+data_num+'.h5','r')
        data_dirc={
                    'den':'den','ene':'ene',
                    'xc':'xc1','yc':'xc2','zc':'xc3',
                    'xl':'xl1','yl':'xl2','zl':'xl3',
                    'xr':'xr1','yr':'xr2','zr':'xr3',
                    'dx':'dx1','dy':'dx2','dz':'dx3',
                    'momx':'momx','momy':'momy','momz':'momz',
                    'sgfx':'sgfx','sgfy':'sgfy','sgfz':'sgfz',
                    'bxr':'bxr','byr':'byr','bzr':'bzr',
                    'bxl':'bxl','byl':'byl','bzl':'bzl',

                    'ndim':'ndim','nvar':'nvar','variable':'varialbe','ndom':'nMesh','nbuf':'nbuf',
                    'lbdry':'leftBdry','rbdry':'rightBdry','t':'t',
                    'CFL':'CFL','boundaryType':'boundaryType','coordType':'coordType','solverType':'solverType','limiterType':'limiterType',
                    'GravConst':'GravConst','snd':'snd','adiGamma':'adiGammas'
                  }
        v_dirc={'vx':'momx','vy':'momy','vz':'momz'}

        for keys,value in data_dirc.iteritems():
            try:
                setattr(self,keys,file[value][:])
            except KeyError:
                setattr(self,keys,None) 
        for keys,values in v_dirc.iteritems():
            try:
                setattr(self,keys,getattr(self,values)/self.den)
            except NameError:
                setattr(self,keys,None)
        self.nmeshx=int(self.ndom[0]+2*self.nbuf)
        self.nmeshy=int(self.ndom[1]+2*self.nbuf)
        if self.ndim ==3:
            self.X,self.Y,self.Z=np.meshgrid(self.xc, self.yc, self.zc)
        else:
            self.X,self.Y=np.meshgrid(self.xc, self.yc)
        file.close()
    def plt_contourf(self,ax,Z,N,cont_argv,*label):
        X=self.X/l_unit
        Y=self.Y/l_unit
        CS = ax.contourf(X,Y,Z,N,**cont_argv)
        try:
            ax.set_xlabel(label[0])
            ax.set_ylabel(label[1])
        except:
            pass
        ax.tick_params(direction='out', length=3, width=0.6, labelsize=6)
        axins=colorbar_top(ax)
        try:
            cb_max=np.max(cont_argv['levels'])
            cb_min=np.min(cont_argv['levels'])
            cb=plt.colorbar(CS,axins,orientation='horizontal',ticks=np.arange(cb_min,cb_max+1,1))
        except KeyError:
            cb=plt.colorbar(CS,axins,orientation='horizontal')
        cb.ax.xaxis.set_ticks_position('top')
        cb.ax.tick_params(direction='out', length=2, width=0.4,labelsize=5)

    #def plt_quiver(ax,X,Y,)

    def plt_b(self,ax):
        nmeshx=self.nmeshx
        nmeshy=self.nmeshy
        dx=self.dx
        dy=self.dy
        bxl=self.bxl
        bxr=self.bxr
        byl=self.byl
        byr=self.byr
        xl=self.xl
        xr=self.xr
        yl=self.yl
        yr=self.yr
        Az = np.zeros((nmeshx+1,nmeshy+1),dtype=float)
        for i in range(1,nmeshx+1):
            Az[i,0]=Az[i-1,0]-byl[0,i-1]*dx[i-1]
        
        for j in range(1,self.nmeshy+1):
            for i in range(0,nmeshx):
                Az[i,j]=Az[i,j-1]+bxl[j-1,i]*dy[j-1]
        
        for j in range(1,nmeshy+1):
            Az[nmeshx,j]=Az[nmeshx,j-1]+bxr[j-1,nmeshx-1]*dy[j-1]
        
        xl = np.append(xl,xr[-1])
        yl = np.append(yl,yr[-1])
        Xl,Yl = np.meshgrid(xl/l_unit,yl/l_unit)
        Q = ax.contour(Xl,Yl,Az.conj().T,5,colors='w',linewidths=0.6,linestyles='solid',levels=np.linspace(np.min(Az),np.max(Az),10))
    def set_range(self,ax,lim):
        ax.set_xlim(-lim,lim)
        ax.set_ylim(-lim,lim)

    def set_mask(self,i):
        dx=self.dx
        dy=self.dy
        dr=np.hypot(dx,dy)
        R=np.hypot(self.X,self.Y)
        r_up=np.hypot(self.xc[i],self.yc[i])+dr/2
        r_low=np.hypot(self.xc[i],self.yc[i])-dr/2
        mask1 = R <= r_up
        mask2 = R >= r_low
        mask = np.logical_and(mask1, mask2)
        return mask
    def cir_avr(self,data):
        cir_avr_data=[]
        for i in range(0,self.nmeshx):
            cir_avr_data.append(np.mean(data[self.set_mask(i)]))
        return np.array(cir_avr_data)
    def check_div_b(self,ax):
        div_b =-self.bxl+self.bxr-self.byl+self.byr
        CS = ax.contourf(self.X, self.Y, div_b)
        plt.colorbar(CS)
                    
    def calc_v(self):
        theta=np.arctan2(self.Y,self.X)
        vr = np.cos(theta)*self.vx+np.sin(theta)*self.vy
        vt = -np.sin(theta)*self.vx+np.cos(theta)*self.vy
        return vr,vt

    def calc_ene(self):
        try:
            vsq=self.vx**2.+self.vy**2.+self.vz**2.
            bx=0.5*(self.bxr+self.bxl)
            by=0.5*(self.byr+self.byl)
            bz=0.5*(self.bzr+self.bzl)
            bsq=bx**2.+by**2.+bz**2.
        except NameError:
            vsq=self.vx**2.+self.vy**2.
            bx=0.5*(self.bxr+self.bxl)
            by=0.5*(self.byr+self.byl)
            bsq=bx**2.+by**2.
        Ek=0.5*self.den*vsq
        Eb=0.5*bsq
        return Ek,Eb
