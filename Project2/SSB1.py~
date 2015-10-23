# importing useful libraries
import numpy as np # numerical package
import matplotlib.pyplot as plt # plotting package
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing


# Definition of constants and variables
keV	= 8.61734E-5 	#Boltzmann constant in eV/K
kerg	= 1.380658E-16	#Boltzmann constant in erg K 
kJoule  = kerg*1e-7	#Boltzmann constant in joule/K
h	= 6.62607E-27	#Planck constant in erg s
c	= 2.99792E10	#Speed of light in cm/s
m_e	= 9.109390E-28	#Electron mass in grams
m_H     = 1.67352E-24   #Hydrogen mass in grams
m_He    = 3.98*m_H 	#Helium mass in grams
m_p     = 1.67262E-24	#Proton mass in grams


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
"""
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
#plt.xscale('log')
#plt.show()
"""
#Find average c value
g_S  = np.average(ptot/colm)
print 'g_S = ',g_S

#Plot ratio hydrogen mass density/total mass density vs height
rho_H  = nhyd * m_H	#Hydrogen mass density
nhel   = nhyd*.1	#Helium number density
rho_He = nhel*m_He	#Helum mass density
"""
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
"""
#Metal fraction
metalfraction = 1. - np.average((rho_H+rho_He)/dens)
print 'Metalfraction = ',metalfraction
"""
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
"""
#Scale height
rhodive = dens[-1]/np.exp(1)
#rhodive = dens[np.where(h==0)][0]/np.exp(1)
rhoe = np.zeros(np.size(h))
rhoe.fill(rhodive)
#Plot gas density vs height
plt.figure(7)
plt.plot(h,dens)
plt.plot(h,rhoe)
#plt.legend([])
plt.title(r'Gas density vs height')
#plt.ylim([0,1.1])
#plt.yscale('log')
plt.ylabel(r'Gas density $\rho$ [g cm$^{-3}$]',size=14)
plt.xlabel(r'Height [km]',size=14)
plt.legend(['Gas density',r'$\rho/e$'])
plt.show()

#Compute gas pressure and plot against height
pgas = pgasptot*ptot #Gas pressure
idealgas = (nhyd + nel)*kerg*temp
"""
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
"""

#Add n_helium and plot again
idealgas_hel = (nhyd + nel+ nhel)*kerg*temp
"""
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
"""
rho_e = nel*m_e	#elecron density
rho_p = nprot*m_p #proton density
"""
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
"""
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
"""
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
"""
#Plot pressure and density stratifications
#together in normalized units in one graph

edens = 10.**elogdens #convert from log to actual
eN = 10.**elogN	#Convert from log to actual
"""
plt.figure(18)
plt.plot(eh,edens/np.max(edens)) #This probably needs a fix
plt.plot(eh,eN/np.max(eN))
plt.yscale('log') #plot on logscale
plt.legend([r'$\rho/\rho_{max}$','$N/N_{max}$'])
plt.title(r'Density and pressure stratification')
plt.ylabel(r'x')
plt.xlabel(r'Height [km]')
#plt.show()
"""
#Plot mean molecular weight vs height
my_E = 10.**(elogdens)/(10.**(elogN)*m_H) #mean molecular weight

"""
plt.figure(19)
plt.plot(eh,my_E)
plt.title(r'Mean molecular weight vs height')
plt.ylabel(r'$\mu_E$')
plt.xlabel(r'Height [km]')
plt.show()
"""
#Estimate density scale height of lower terrestrial atmosphere
eH_p = -eh[2] / np.log(edens[2]/edens[np.where(eh==0)][0])
print 'Scale height H_p low earth atmosphere =  ',eH_p,' km'

#Compare terrestrial parameter values to solar
print '------------------'
print 'All values calculated at h = 0'
print 'Matter density earth = ',dens[np.where(h==0)][0]
print 'Matter density sun   = ',edens[np.where(eh==0)][0]

N = (nhyd + nhel + nel)
print 'Particle density sun   = ',N[np.where(h==0)][0]
print 'Particle density earth = ',eN[np.where(eh==0)][0]

eptot = 10**elogP

print 'Pressure sun   = ', ptot[np.where(h==0)][0]
print 'Pressure earth = ', eptot[np.where(eh==0)][0]

print 'Temperature sun   = ', temp[np.where(h==0)][0]
print 'Temperature earth = ', etemp[np.where(eh==0)][0]

print 'Ratio particle densities = ', eN[np.where(eh==0)][0]/N[np.where(h==0)][0]
print '------------------'

#Estimate atmospheric column mass earth
g_E = 980.665 #cm s^-2
ecolm = eptot[np.where(eh==0)][0]/g_E
print 'Column mass earth surface = ', ecolm
print 'Column mass sun surface = ', colm[np.where(h==0)][0]

#Compare N_phot to particle density in the air on earth
#N_phot = 
#print 'Ratio photons and particles earth atmosphere = ',N_phot/eN


#Part 2 Observed solar continua

#Write code to read table 5 (solspect.dat)
solspect = np.loadtxt('solspect.dat',unpack=True)
wavelength = solspect[0,:]

F_lambda = solspect[1,:]
F_lambda_c = solspect[2,:]
I_lambda = solspect[3,:]
I_lambda_c = solspect[4,:]

#Plot the four distributions in one fig
#over the range wavelength = 0-2 micro meter.
"""
plt.figure(20)
plt.plot(wavelength,F_lambda)
plt.plot(wavelength,F_lambda_c)
plt.plot(wavelength,I_lambda)
plt.plot(wavelength,I_lambda_c)
plt.xlim([0,2])
plt.title('The four spectral distributions')
plt.legend([r'$F_\lambda$',r'$F_\lambda^c$',r'$I_\lambda$',r'$I_\lambda^c$'])
plt.xlabel(r'Wavelength $\lambda$ [$\mu m$]')
plt.ylabel(r'Distributions [erg cm$^{-2}$s$^{-1}$ster$^{-1}\mu$m$^{-1}$]')
#plt.show()
"""
#Check that continiuum intensity reaches 4.6*10^10 at max
print 'I_lambda_c = %e'%np.max(I_lambda_c) #check okay

#Convert distributions into values per frequency
F_nu 	= F_lambda 	/c* (wavelength*1e-4)**2 *1e4 
F_nu_c 	= F_lambda_c 	/c* (wavelength*1e-4)**2 *1e4
I_nu 	= I_lambda 	/c* (wavelength*1e-4)**2 *1e4
I_nu_c 	= I_lambda_c 	/c* (wavelength*1e-4)**2 *1e4

frequency = c/(wavelength*1e-4)
"""
plt.figure(21)
plt.plot(wavelength,F_nu)
plt.plot(wavelength,F_nu_c)
plt.plot(wavelength,I_nu)
plt.plot(wavelength,I_nu_c)
plt.title('The four spectral distributions')
plt.legend([r'$F_\nu$',r'$F_\nu^c$',r'$I_\nu$',r'$I_\nu^c$'])
plt.xlabel(r'Wavelength $\lambda$ [$\mu m$]')
plt.ylabel(r'Distribtuions [erg cm$^{-2}$s$^{-1}$ster$^{-1}$Hz$^{-1}$]')
#plt.show()
"""
print 'Peak I_nu_c = %e'%np.max(I_nu_c)

#Use planck func form SSA
def planck(temp,wav):
	'''This function takes in wavelength in micrometers and temperature in K.
	and returns the planck function with units [1e10 erg/cm^2/s/micrometer/steradian]'''

	h	= 6.62607E-27	#Planck constant in erg s
	B = 2.*h*c*c/((wav*1e-4)**5) / (np.exp( (h*c) / ((wav*1e-4)*kerg*temp) ) - 1 )
	return B*1e-14 #1e-14 coverts to same units as in tables
"""
planckfit = planck(6450,wavelength) 
plt.figure(22)
plt.plot(wavelength,I_lambda_c)
plt.plot(wavelength,planckfit)
plt.title(r'Observed intensity of continuum vs Planck fit')
plt.ylabel(r'Intensity [$10^{10}$ erg cm$^{-2}$s$^{-1}$ster$^{-1}\mu$m$^{-1}$]')
plt.xlabel(r'Wavelength $\lambda$ [$\mu m$]')
plt.legend([r'$I_\lambda^c$',r'$B_\lambda(6450)$'])
plt.show()
"""
