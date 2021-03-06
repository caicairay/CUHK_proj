#imports
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import hdf_class as hc
import pickle
import sys

save_fig=False
switcher=1

data_num=sys.argv[1]
DATA=hc.HDF5_DATA(data_num=data_num)

#f=open("../init_data/init_ek_log.txt","r")
#Ek_init=pickle.load(f)
#f.close()
#f=open("../init_data/init_eb_log.txt","r")
#Eb_init=pickle.load(f)
#f.close()
#f=open("../init_data/init_vt_log.txt","r")
#vt_init=pickle.load(f)
#f.close()

if save_fig==True:
    switcher=int(sys.argv[2])

def set_mask(lim):
    R=np.hypot(DATA.X,DATA.Y)
    mask = R/hc.l_unit <= lim
    return mask
#density & B-field#
if switcher == 1:  
    #__Control Pannel__#
    lim=200
    cont_argv={'cmap':hc.custom_cm(),
               'levels':np.linspace(-11,-8,100),
               'extend':'both' 
              }
    
    #__Figure Setup__#
    fig,ax=plt.subplots(1,1)
    fig.suptitle("Time: "+str(DATA.t/hc.t_unit)+"yr")
    DATA.plt_contourf(ax,np.log10(DATA.den),0,cont_argv,'AU')
    DATA.plt_b(ax)
    DATA.set_range(ax,lim)
    #ax.set_xlim(-40,-80)  #for zoom in only
    #ax.set_ylim(110,140) #
    #ax.grid(True,  linestyle='--', color='r') #which='minor', axis='both',
    ax.set_aspect('equal')

#Ek & Eb#
elif switcher == 2:
    #__Control Pannel__#
    lim=200
    cont_argv1={'cmap':cm.coolwarm,
               'levels':np.linspace(0,2,20),
               'extend':'both' 
              }
    cont_argv2={'cmap':cm.coolwarm,
               'levels':np.linspace(-4,-2,20),
               'extend':'both' 
              }
    
    #__Figure Setup__#
    Ek,Eb=DATA.calc_ene()
    fig,ax=plt.subplots(1,2)
    fig.suptitle("Time: "+str(DATA.t/hc.t_unit)+"yr")
    DATA.plt_contourf(ax[0],np.log10(Ek),30,cont_argv1,'AU')
    DATA.plt_contourf(ax[1],np.log10(Eb),30,cont_argv2,'AU')
    DATA.plt_b(ax[0])
    DATA.plt_b(ax[1])
    DATA.set_range(ax[0],lim)
    DATA.set_range(ax[1],lim)
    ax[0].set_aspect('equal')
    ax[1].set_aspect('equal')

#Ek Eb avr#
elif switcher == 3:
    #__Control Pannel__#
    lim=200
    
    #__Figure Setup__#
    R=np.hypot(DATA.xc,DATA.yc)
    Ek,Eb=DATA.calc_ene()
    Ek_avr=DATA.cir_avr(np.log10(Ek))
    Eb_avr=DATA.cir_avr(np.log10(Eb))
    Ekl,Ekr=np.split(Ek_avr,2) 
    Ebl,Ebr=np.split(Eb_avr,2) 
    Rl,Rr=np.split(R/hc.l_unit,2)
    fig,ax=plt.subplots(1,1)
    fig.suptitle("Time: "+str(DATA.t/hc.t_unit)+"yr")
    ax.plot(Rr,Ekr,'r',label="Ek")
    ax.plot(Rr,Ek_init,'--k')
    ax.plot(Rr,Ebr,'y',label="Eb")
    ax.plot(Rr,Eb_init,'--k')
    ax.set_ylim(-5,2)
    ax.set_xlim(0,lim)
    ax.set_yticks(np.arange(-5, 3, step=1))
    ax.legend(loc=1)

#vt avr#
elif switcher == 4:
    #__Control Pannel__#
    lim=200
    
    #__Figure Setup__#
    R=np.hypot(DATA.xc,DATA.yc)
    vr,vt=DATA.calc_v()
    vr_avr=DATA.cir_avr(np.log10(vr))
    vt_avr=DATA.cir_avr(np.log10(vt))
    vrl,vrr=np.split(vr_avr,2) 
    vtl,vtr=np.split(vt_avr,2) 
    Rl,Rr=np.split(R/hc.l_unit,2)
    fig,ax=plt.subplots(1,1)
    fig.suptitle("Time: "+str(DATA.t/hc.t_unit)+"yr")
    ax.plot(Rr,vtr,'y',label="log10(vt)")
    ax.plot(Rr,vt_init,'--k')
    #ax.plot(Rr,vrr,'r',label="vr(km/s)")
    #ax.plot(Rr,vr_init,'--k')
    ax.set_ylim(5,6)
    ax.set_xlim(0,lim)
    ax.set_yticks(np.arange(5, 7, step=.5))
    ax.legend(loc=1)

elif switcher == 5:
    #__Control Pannel__#
    lim=200
    cont_argv1={'cmap':cm.RdBu_r,
               #'levels':np.linspace(-4,4,20),
               #'extend':'both' 
              }
    
    #__Figure Setup__#
    E = DATA.ene
    Ek,Eb=DATA.calc_ene()
    Eu = E-Ek-Eb
    #fig,ax=plt.subplots(1,1)
    #fig.suptitle("Time: "+str(DATA.t/hc.t_unit)+"yr")
    #DATA.plt_contourf(ax,E,30,cont_argv1,'AU')
    #DATA.plt_b(ax)
    #DATA.set_range(ax,lim)
    #ax.set_aspect('equal')
    print np.sum(Ek[set_mask(10)])

#_______________________#
if save_fig==True:
    fig_path = sys.argv[3]
    plt.savefig(fig_path+"/hdf_"+str('{:0>4}'.format(data_num))+".png",dpi=250)
else:
    print 'saving'
    plt.savefig("./test"+str(switcher)+".png",dpi=200)
