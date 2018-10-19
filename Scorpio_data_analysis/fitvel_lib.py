### imports ###
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
import sys

plt_case = 1  # 1 for density , 2 for integrated mass
enable_vel_mom = 1 # 0 for close, 1 for velosity, 2 for momentum
### Constants ###
G = 4.3e-3 #pc/Msun (km/s)^2
pc_to_AU = 206.2634 #1pc = 206.2634 kAU
#G = G * pc_to_AU #kAU/Mun (km/s)^2 
G = G / pc_to_AU #pc^2/Mun/kAU (km/s)^2 

### Set upperbdry ###
upperbdry = np.linspace(0.001, 4.0, 1000) # kAU

### Set disk density profile ###
def den_profile(R):
    den = 0.1*np.exp(np.log(6./5.)/(-1.0)*R**2.0)  #M_sun/pc^2
    return den

### Calculate total disk mass ###
def disk_mass(upperbdry):  # [M_sun/pc^2*kAU^2]
    dmass = []
    for i in upperbdry:
        dmass_ = integrate.quad(lambda R: 2*np.pi*R*den_profile(R), 0 ,i)
        dmass.append(dmass_[0])
    print max(dmass)
    return dmass

### Calculate velocity profile ###
def vel_need(upperbdry): # [km/s] 
    v = []
    dmass = disk_mass(upperbdry)
    for i in range(0, len(upperbdry)):
        v_ = np.sqrt(G*dmass[i]/upperbdry[i])
        v.append(v_)
    return v
def vel_uplim(upperbdry): 
    v = []
    dmass = disk_mass(upperbdry)
    for i in range(0, len(upperbdry)):
        v_ = np.sqrt(2*G*dmass[i]/upperbdry[i])
        v.append(v_)
    return v

### Fitting velocity profile ###
def fit_vel(upperbdry):
    x = upperbdry
    y = []
    for i,j in zip(vel_need(upperbdry),vel_uplim(upperbdry)):
        y_ = (i+j)/2
        y.append(y_)
    popt,pcov = curve_fit(lambda r, a,b: a*np.sqrt(r/(b+r)),x,y)
    vel_profile = popt[0]*np.sqrt(upperbdry/(popt[1]+upperbdry))
    print popt
    return vel_profile

### momentum ###
def cal_mom(upperbdry):
    den = []
    for i in range(0,len(upperbdry)):
        den_=den_profile(upperbdry[i])
        den.append(den_)
    return fit_vel(upperbdry)*den

### ploting ###
def plt_curve():
    fig, ax1 = plt.subplots()
    if plt_case==1:
        ax1.plot(upperbdry, den_profile(upperbdry), 'b')
    if plt_case==2:
        ax1.plot(upperbdry, disk_mass(upperbdry),'g')
    ax1.set_xlabel('R')
    ax1.set_ylabel('den', color='b')
    ax1.tick_params('y', colors='b')
    if enable_vel_mom == True:
        ax2 = ax1.twinx()
        if enable_vel_mom ==1:
            ax2.plot(upperbdry, vel_need(upperbdry), 'r')
            ax2.plot(upperbdry, vel_uplim(upperbdry), 'r--')
            ax2.plot(upperbdry, fit_vel(upperbdry), 'g--')
            ax2.set_ylabel('Velocity', color='r')
        if enable_vel_mom ==2:
            ax2.plot(upperbdry, cal_mom(upperbdry), 'r')
            ax2.set_ylabel('momentum',color='r')
        ax2.tick_params('y', colors='r')
    fig.tight_layout()
    plt.show()

if __name__=="__main__":
    plt_curve()
