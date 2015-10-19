# importing useful libraries (you may need more)
import numpy as np # numerical package
import matplotlib.pyplot as plt # plotting package
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing


# Definition of constants and variables
keV	= 8.61734E-5 	#Boltzmann constant in eV/deg
kerg	= 1.380658E-16	#Boltzmann constant in erg K 
h	= 6.62607E-27	#Planck constant in erg s
c	= 2.99792E10	#Speed of light in cm/s
m_e	= 9.109390E-28	#Electron mass in grams
m_H     = 1.67352E-24   #Hydrogen mass in grams
m_He    = 3.98*m_H 	#Helium mass in grams
m_p     = 1.67262E-24	#Proton mass in grams

#Define planck function
def planck(temp,wav):
	B = 2.*h*c*c/(wav**5) / (np.exp( (h*c) / (wav*kerg*temp) ) - 1 )
	return B
#Part 1
#Read in falc data and arrange in arrays
falc = np.loadtxt('falc.dat',unpack=True)
h        = np.array(falc[0,:])
tau5     = np.array(falc[1,:])
colm     = np.array(falc[2,:])
temp     = np.array(falc[3,:])
vturb    = np.array(falc[4,:])
nhyd     = np.array(falc[5,:])
nprot    = np.array(falc[6,:])
nel      = np.array(falc[7,:])
ptot     = np.array(falc[8,:])
pgasptot = np.array(falc[9,:])
dens     = np.array(falc[10,:])
#Plot falc model
plt.figure(0)
plt.plot(h,temp)
plt.ylim([3000,10000])
plt.title(r'FALC Model')
plt.xlabel(r'Height [km]',size=14)
plt.ylabel(r'Temperature [K]',size=14)

#1.2
#Plot total pressure vs column mass linearly
plt.figure(1)
plt.plot(colm,ptot)
plt.title(r'$P_{total}$ vs column mass m')
plt.ylabel(r'$P_{total}$ [dyn cm$^{-2}$]',size=14)
plt.xlabel(r'Column mass m [$g cm^{-2}$]',size=14)
#logarithimically
plt.figure(2)
plt.plot(colm,ptot)
plt.title(r'$P_{total}$ vs column mass m logarithmically')
plt.ylabel(r'$P_{total}$ [dyn cm$^{-2}$]',size=14)
plt.xlabel(r'Column mass m [g cm$^{-2}$]',size=14)
plt.yscale('log')
plt.xscale('log')
#plt.show()

#Find average c value
c  = ptot/colm
print 'c = ',np.average(c)

#Plot ratio hydrogen mass density/total mass density vs height
rho_H  = nhyd * m_H	#Hydrogen mass density
nhel   = nhyd*.1	#Helium number density
rho_He = nhel*m_He	#Helum mass density

plt.figure(4)
plt.plot(h,rho_H/dens)
plt.plot(h,rho_He/dens)
plt.plot(h,(rho_H+rho_He)/dens)

plt.legend([r'$\rho_{H}/\rho_{total}$',r'$\rho_{He}/\rho_{total}$',r'$(\rho_{H}+\rho_{He})/\rho_{total}$'],loc='best')
plt.title(r'Density fractions of total as a function of height')
plt.ylim([0,1.1])
plt.ylabel(r'$\rho/\rho_{total}$',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()

#Metal fraction
metalfraction = 1. - np.average((rho_H+rho_He)/dens)
print 'Metalfraction = ',metalfraction

#Plot column mass vs height
plt.figure(5)
plt.plot(h,colm)
#plt.legend([])
plt.title(r'Column mass vs height')
#plt.ylim([0,1.1])
plt.ylabel(r'Column mass [g cm$^{-2}$]',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()

#Logarithmic
plt.figure(6)
plt.plot(h,colm)
#plt.legend([])
plt.title(r'Column mass vs height logarithmic')
#plt.ylim([0,1.1])
plt.yscale('log')
plt.ylabel(r'Column mass [g cm$^{-2}$]',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()

#Plot gas density vs height
plt.figure(7)
plt.plot(h,dens)
#plt.legend([])
plt.title(r'Gas density vs height')
#plt.ylim([0,1.1])
#plt.yscale('log')
plt.ylabel(r'Gas density $\rho$ [g cm$^{-3}$]',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()


#Compute gas pressure and plot against height
pgas = pgasptot*ptot #Gas pressure
idealgas = (nhyd + nel)*kerg*temp
#Overplot product (nhyd + nel)kT
plt.figure(8)
plt.plot(h,pgas) #from falc
plt.plot(h,idealgas) #from idealgas law
plt.legend(['Falc','Ideal gas'])
plt.title(r'Gas pressure vs height')
#plt.ylim([0,1.1])
#plt.yscale('log')
plt.ylabel(r'Gas pressure $P_{gas}$ [dyn cm$^{-2}$]',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()

#Plot ratio between curves
plt.figure(9)
plt.plot(h,pgas/idealgas) #ratio curves
plt.title(r'Ratio gas pressure from falc and ideal gas law')
plt.ylabel(r'Gas pressure $P_{gas}$ [dyn cm$^{-2}$]',size=14)
plt.xlabel(r'Height [km]',size=14)
plt.show()


#Add n_helium and plot again
idealgas_hel = (nhyd + nel+ nhel)*kerg*temp
#Overplot product (nhyd + nel + nhel)kT
plt.figure(10)
plt.plot(h,pgas) #from falc
plt.plot(h,idealgas_hel) #from idealgas with helium
plt.legend(['Falc','Ideal gas with helium'])
plt.title(r'Gas pressure vs height')
#plt.ylim([0,1.1])
#plt.yscale('log')
plt.ylabel(r'$P_{gas}$ [dyn cm$^{-2}$]',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()

#Plot ratio between curves
plt.figure(11)
plt.plot(h,pgas/idealgas_hel) #ratio curves
plt.title(r'Ratio gas pressure from falc and ideal gas law')
plt.ylabel(r'Ratio $P_{FALC}/P_{IG}$',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()

#Plot total hydrogen density vs height
#and overplot electron density,
#proton density, and density of electrons not from
#hydrogen ionization.

rho_e = nel*m_e	#elecron density
rho_p = nprot*m_p #proton density

plt.figure(12)
plt.plot(h,rho_H) #hydrogen vs h
plt.plot(h,rho_e) #electron vs h
plt.plot(h,rho_p) #proton vs h
#plt.plot(h,rho_nie) #ratio curves
#plt.legend([r'\rho_H',r'\rho_e',r'\rho_p',r'\rho_{nie}'])
plt.legend([r'$\rho_H$',r'$\rho_e$',r'$\rho_p$'])
plt.title(r'Select densities vs height')
plt.ylabel(r'\rho [g cm$^{-3}$]',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()

#Plot ioniaztion fraction of hydrogen
#logarithimically against height
"""
hyd_ion = 

plt.figure(13)
plt.plot(hyd_ion,h)
plt.title(r'Ioniaztion fraction of hydrogen vs height')
plt.ylabel(r'Ionization fraction',size=14)
plt.xlabel(r'Height [km]',size=14)
plt.show()
"""

#Write code to read earth.dat
earth = np.loadtxt('earth.dat',unpack=True)
eh = np.array(earth[0,:])
elogP = np.array(earth[1,:])
etemp = np.array(earth[2,:])
elogdens = np.array(earth[3,:])
elogN = np.array(earth[4,:])

#Plot temp,pressure,density,gas density
plt.figure(14)
plt.plot(eh,elogP) #ratio curves
plt.title(r'Air pressure vs height logarithmic')
plt.ylabel(r'log$P [dyn cm$^{-2}$]$',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()

plt.figure(15)
plt.plot(eh,etemp) #ratio curves
plt.title(r'Temperature vs height')
plt.ylabel(r'Temperature [K]',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()

plt.figure(16)
plt.plot(eh,elogdens) #ratio curves
plt.title(r'Gas density vs height logarithmic')
plt.ylabel(r'log $\rho$ [g cm$^{-3}$]',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()

plt.figure(17)
plt.plot(eh,elogN) #ratio curves
plt.title(r'Particle density vs height')
plt.ylabel(r'log $N$ [cm$^{-3}$]',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()

#Plot pressure and density stratifications
#together in normalized units in one graph
plt.figure(18)
plt.plot(eh,elogdens/np.abs(np.min(elogdens))) #This probably needs a fix
plt.plot(eh,elogN/np.max(elogN))
plt.legend([r'$\rho/\rho_{max}$','$N/N_{max}$'])
plt.title(r'Density and pressure stratification')
plt.ylabel(r'x')
plt.xlabel(r'Height [km]')
#plt.show()

#Plot mean molecular weight vs height
my_E = 10.**(elogdens)/(10**(elogN)*m_H) #mean molecular weight
plt.figure(19)
plt.plot(eh,my_E)
plt.title(r'Mean molecular weight vs height')
plt.ylabel(r'$\mu_E$')
plt.xlabel(r'Height [km]')
plt.show()

