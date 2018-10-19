# imports
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt

# Constants
#G = 4*(np.pi**2) # Gravitational constant: AU^3 yr^-2 M_sun^-1
G = 4.3 #pc^2/Msun/kpc (km/s)^2
# pc to AU
snd = 25.
# Set upperbdry
upperbdry = np.linspace(0.001, 15.0, 1000) # kpc

# Set disk density profile
def den_profile(R): #Msun/pc^2
    den = 110.*np.exp(np.log(6.0/5.0)/(-15.0)*R**2)  # density  M_sun/pc^2
    return den

# Calculate total disk mass 
def disk_mass(upperbdry):
    dmass = []
    for i in upperbdry:
        dmass_ = integrate.quad(lambda R: 2*np.pi*R*den_profile(R), 0 ,i)
        dmass.append(dmass_[0])
    print max(dmass)
    return dmass

# Calculate velocity profile, km/s
def vel_need(upperbdry): 
    v = []
    dmass = disk_mass(upperbdry)
    for i in range(0, len(upperbdry)):
        v_ = np.sqrt(G*dmass[i]/upperbdry[i]+2.0*snd**2.0*upperbdry[i]**2.0*np.log(6./5.)/(-15.))
        v.append(v_)
    return v
def vel_uplim(upperbdry): 
    v = []
    dmass = disk_mass(upperbdry)
    for i in range(0, len(upperbdry)):
        v_ = np.sqrt(2*G*dmass[i]/upperbdry[i]+2.0*snd**2.0*upperbdry[i]**2.0*np.log(6./5.)/(-15.))
        v.append(v_)
    return v

# Fitting velocity profile
def fit_vel(upperbdry):
    x = upperbdry
    #y = []
    #for i,j in zip(vel_need(upperbdry),vel_uplim(upperbdry)):
    #    y_ = (i+j)/2
    #    y.append(y_)
    y = vel_uplim(upperbdry)
    popt,pcov = curve_fit(lambda r, a,b: a*np.sqrt(r/(b+r)),x,y)
    #vel_profile = popt[0]*np.sqrt(upperbdry/(popt[1]+upperbdry))
    vel_profile = 220.*np.sqrt(upperbdry/(4.+upperbdry))
    print popt
    return vel_profile
#momentum
def cal_mom(upperbdry):
    den = []
    for i in range(0,len(upperbdry)):
        den_=den_profile(upperbdry[i])
        den.append(den_)
    return fit_vel(upperbdry)*den
#print fit_vel(upperbdry)

def plt_curve():
    #Ploting
    fig, ax1 = plt.subplots()
    ax1.plot(upperbdry, den_profile(upperbdry), 'b')
    ax1.set_xlabel('R')
    ax1.set_ylabel('den', color='b')
    ax1.tick_params('y', colors='b')
    #ax1.plot(upperbdry, disk_mass(upperbdry),'g')
    ax2 = ax1.twinx()
    ax2.plot(upperbdry, vel_need(upperbdry), 'r')
    ax2.plot(upperbdry, vel_uplim(upperbdry), 'r--')
    ax2.plot(upperbdry, fit_vel(upperbdry), 'g--')
    ax2.set_ylabel('Velocity', color='r')
    #ax2.plot(upperbdry, cal_mom(upperbdry), 'r')
    #ax2.set_ylabel('momentum',color='r')
    ax2.tick_params('y', colors='r')
    #plt.ylim(0, 1000)
    
    fig.tight_layout()
    plt.show()
#plt.savefig('/home/zhuo/Desktop/fitted_vel.png',format = 'png', dpi = 300)

if __name__=="__main__":
    plt_curve()
