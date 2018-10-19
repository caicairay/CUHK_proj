import numpy as np
import matplotlib.pyplot as plt
import fitvel_lib as ci
r=np.linspace(0.001,3000,1000); #[AU]
dr=r[2]-r[1]
snd = 25.0 #km/s

#rho = 110.0*np.exp(np.log(6.0/5.0)/(-15.0)*r**2) #[Msun/pc^2]
rho = []
for i in r:
    rho_ = ci.den_profile(i)
    rho.append(rho_)
#rho = ci.den_profile(r) #Msun/AU^2
#print type(rho)
#cv=220.0*np.sqrt(r/(0.4+r));
cv = ci.vel_need(r) #km/s
G=4.3
# pc to AU
pc_to_AU = 206263.4 #1pc = 206263.4 AU
G = G * pc_to_AU #AU/Msun (km/s)^2

omega=cv/r

kappa=np.sqrt(4*omega**2*(1+r/(2.0*omega)*np.gradient(omega,dr)))
Q=[]
for i in rho:
    Q_ = snd*kappa/(np.pi*G*i)
    Q.append(Q_)
#Q = snd*kappa/(np.pi*G*rho)
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
#plt.ylim(0,5)
plt.show()
