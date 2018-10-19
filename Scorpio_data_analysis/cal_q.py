import numpy as np
import matplotlib.pyplot as plt
import fitvel_lib as fl 
r=np.linspace(0.001,8.0,1000); #[kAU]
dr=r[2]-r[1]
snd = 0.6 #km/s

rho = fl.den_profile(r)
cv = fl.fit_vel(r) #km/s
G=4.3e-3
# pc to AU
pc_to_AU = 206.2634 #1pc = 206.2634 kAU
G = G * pc_to_AU #kAU/Msun (km/s)^2

omega=cv/r

kappa=np.sqrt(4*omega**2*(1+r/(2.0*omega)*np.gradient(omega,dr)))
Q = snd*kappa/(np.pi*G*rho)
plt.figure(1)
plt.plot(r,rho)
plt.xlabel("r")
plt.ylabel("rho")

plt.figure(2)
plt.plot(r,cv)
plt.xlabel("r")
plt.ylabel("cv")

plt.figure(3)
plt.subplot(311)
plt.plot(r, omega)
plt.xlabel("r")
plt.ylabel("omega")
plt.subplot(312)
plt.plot(r, omega-kappa/2)
plt.xlabel("r")
plt.ylabel("omega-kappa/2")
plt.subplot(313)
plt.plot(r, omega+kappa/2)
plt.xlabel("r")
plt.ylabel("omega+kappa/2")

plt.figure(4)
plt.plot(r, Q)
plt.xlabel("r")
plt.ylabel("Q")
plt.ylim(0,100)
plt.show()
