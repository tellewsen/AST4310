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
m_u	= 1.660538e-24	#Atomic mass unit in grams
R_sun   = 696300	#Radius sun in km
D_sun   = 1.496e8 	#Distance sun-earth in km
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
#1.1
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
"""
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
"""
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
#rho_e = nel*m_e	#electron density
#rho_p = nprot*m_p #proton density
n_eb = (nhyd - nprot)
"""
plt.figure(12)
plt.plot(h,nhyd) #hydrogen vs h
plt.plot(h,nel) #electron vs h
plt.plot(h,nprot) #proton vs h
plt.plot(h,n_eb) #bound electron vs h
plt.legend([r'$n_H$',r'$n_e$',r'$n_p$',r'$n_{be}$'])
plt.title(r'Select densities vs height')
plt.ylabel(r'Number density [cm$^{-3}$]',size=14)
plt.xlabel(r'Height [km]',size=14)
plt.yscale('log')
plt.show()
"""
#Plot ioniaztion fraction of hydrogen
#logarithimically against height

hyd_ion = nprot/nhyd
"""
plt.figure(13)
plt.plot(h,hyd_ion)
plt.title(r'Ionization fraction of hydrogen vs height')
plt.ylabel(r'Ionization fraction',size=14)
plt.xlabel(r'Height [km]',size=14)
plt.yscale('log')
plt.ylim([0,1.1])
plt.xlim([-100,2220])
plt.show()
"""
#Photon density
N_phot = 20*temp**3

print 'N_phot for deep photosphere = ',N_phot[np.where(h==np.min(h))][0]
print 'N_H for deep photosphere    = ',nhyd[np.where(h==np.min(h))][0]
Teff = 5770.
N_phot_high = 20.*Teff**3/(2.*np.pi)
print 'N_phot for highest point sun = %e'%N_phot_high
print 'N_H for highest point sun    = %e'%nhyd[np.where(h==np.max(h))][0]

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
plt.show()

plt.figure(17)
plt.plot(eh,elogN) #ratio curves
plt.title(r'Particle density vs height')
plt.ylabel(r'log $N$ [cm$^{-3}$]',size=14)
plt.xlabel(r'Height [km]',size=14)
#plt.show()
"""
#Plot pressure and density stratifications
#together in normalized units in one graph

eN = 10.**elogN	#Convert from log to actual
edens = 10.**elogdens #convert from log to actual

edense = np.zeros(np.size(eh))
edense.fill(edens[np.where(eh==0)][0]/np.exp(1))
"""
plt.figure(17)
plt.plot(eh,edens)
plt.plot(eh,edense)
plt.show()
"""
"""
plt.figure(18)
plt.plot(eh,edens/np.max(edens))
plt.plot(eh,eN/np.max(eN))
plt.yscale('log') #plot on logscale
plt.legend([r'$\rho/\rho_{max}$','$N/N_{max}$'])
plt.title(r'Density and pressure stratification')
plt.ylabel(r'x')
plt.xlabel(r'Height [km]')
#plt.show()
"""
#Plot mean molecular weight vs height
my_E = edens/(eN*m_H) #mean molecular weight

"""
plt.figure(19)
plt.plot(eh,my_E)
plt.title(r'Mean molecular weight vs height')
plt.ylabel(r'$\mu_E$')
plt.xlabel(r'Height [km]')
plt.show()
"""
#Estimate density scale height of lower terrestrial atmosphere
g_E = 980.665 #cm s^-2
eH_p = kJoule*etemp/((my_E*m_u*1e-3)*(g_E*1e-2))*1e-3
print 'Scale height H_p low earth atmosphere =  ',eH_p[0],' km'

#Compare terrestrial parameter values to solar
print '------------------'
print 'All values calculated at h = 0'

print 'Matter density sun   = ',dens[np.where(h==0)][0]
print 'Matter density earth = ',edens[np.where(eh==0)][0]

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
N_phot_sun = np.pi*R_sun**2/D_sun**2*N_phot_high
N_phot_earth = 20*etemp[np.where(eh==0)][0]**3
print 'Photons from sun: %e'%N_phot_sun
print 'Photons from earth:%e'%N_phot_earth
print 'Particle density surface of earth: %e'%eN[np.where(eh==0)][0]
print 'Ratio photons from sun and particles earth atmosphere = ',N_phot_sun/eN[np.where(eh==0)][0]
print 'Ratio photons from sun and photons from earth atmosphere = ',N_phot_sun/N_phot_earth


#Part 2 Observed solar continua

#Write code to read table 5 (solspect.dat)
solspect = np.loadtxt('solspect.dat',unpack=True)
wavelength = solspect[0,:]

F_lambda = solspect[1,:]#*1e10
F_lambda_c = solspect[2,:]#*1e10
I_lambda = solspect[3,:]#*1e10
I_lambda_c = solspect[4,:]#*1e10

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
conversionfactor = 1./c* (wavelength*1e-4)**2 *1e4
F_nu 	= F_lambda 	*conversionfactor
F_nu_c 	= F_lambda_c 	*conversionfactor
I_nu 	= I_lambda 	*conversionfactor
I_nu_c 	= I_lambda_c 	*conversionfactor

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
plt.xlim([0,2])
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

def planck_nu(temp,wav):
	'''This function takes in wavelength in micrometers and temperature in K.
	and returns the planck function with units [1e10 erg/cm^2/s/micrometer/steradian]'''
	h	= 6.62607E-27	#Planck constant in erg s
	B = 2.*h*c*c/((wav*1e-4)**5) / (np.exp( (h*c) / ((wav*1e-4)*kerg*temp) ) - 1 )
	return B*1e-14*conversionfactor #1e-14 coverts to same units as in tables


planckfit1 = planck(6200,wavelength) 
planckfit2 = planck(6300,wavelength) 
planckfit3 = planck(6400,wavelength)
"""
plt.figure(22)
plt.plot(wavelength,I_lambda_c)
plt.plot(wavelength,planckfit1)
plt.plot(wavelength,planckfit2)
plt.plot(wavelength,planckfit3)
plt.xlim([0,2])
plt.title(r'Observed intensity of continuum vs Planck fit')
plt.ylabel(r'Intensity [$10^{10}$ erg cm$^{-2}$s$^{-1}$ster$^{-1}\mu$m$^{-1}$]')
plt.xlabel(r'Wavelength $\lambda$ [$\mu m$]')
plt.legend([r'$I_\lambda^c$',r'$B_\lambda(6200)$',r'$B_\lambda(6300)$',r'$B_\lambda(6400)$'])
#plt.show()
"""
#planckfit for frequency
planckfit4 = planck_nu(6200,wavelength)*1e15 #1e15 changes units on y axis.
planckfit5 = planck_nu(6300,wavelength)*1e15 #note the change from 10^10
planckfit6 = planck_nu(6400,wavelength)*1e15 #to 10^-5.
"""
plt.figure(23)
plt.plot(wavelength,I_nu_c*1e15)
plt.plot(wavelength,planckfit4)
plt.plot(wavelength,planckfit5)
plt.plot(wavelength,planckfit6)
plt.xlim([0,2])
plt.title(r'Observed intensity of continuum vs Planck fit')
plt.ylabel(r'Intensity [$10^{-5}$ erg cm$^{-2}$s$^{-1}$ster$^{-1}$Hz$^{-1}$]')
plt.xlabel(r'Wavelength $\lambda$ [$\mu m$]')
plt.legend([r'$I_\nu^c$',r'$B_\nu(6200)$',r'$B_\nu(6300)$',r'$B_\nu(6400)$'])
plt.show()
"""
#Invert planckfunction

def planck_invert(I,wavelength):
	'''Takes in intensity in units 10^10 erg cm^-2 s^-1 ster^-1 
	\mu m^-1, and wavelength in units \mu m, and returns brightness 
	temperatre T_b '''

	h	= 6.62607E-27	#Planck constant in erg s
	T_b = h*c/(wavelength*1e-4*kerg) *1./ np.log((2.*h*c*c)/(I*1e14*(wavelength*1e-4)**5) +1.)
	return T_b

T_b  = planck_invert(I_lambda_c,wavelength)#Brightness temperature solar continuum
"""
plt.figure(24)
plt.plot(wavelength,T_b)
plt.title(r'Brightness temperature vs wavelength for the solar continuum')
plt.ylabel(r'Brightness temperature $T_b$[K] ')
plt.xlabel(r'Wavelength [$\mu$m]')
plt.show()
"""

#2.2 Continuous extinction
def exthmin(wav,temp,eldens):
	"""
	in: 	wav = wavelength [Angstrom] (float or fltarr)
		temp = temperature [K]
		eldens = electron density [electrons cm-3]
	
	out: 	H-minus bf+ff extinction [cm^2 per neutral hydrogen atom]
		assuming LTE ionization H/H-min
	"""

	h = 6.626076E-27   #Planck constant (erg s) 
	theta = 5040. / temp
	elpress = eldens*kerg*temp
	sigmabf = ( 1.99654 -1.18267E-5*wav +2.64243E-6*wav**2 
	-4.40524E-10*wav**3+3.23992E-14*wav**4 -1.39568E-18*wav**5 
	+2.78701E-23*wav**6 )

	sigmabf *= 1e-18	
	for i in range(np.size(wav)-1):
		if sigmabf[i] >= 16444:
			sigmabf[i] = 0

	graysaha=4.158E-10*elpress*theta**2.5*10.**(0.754*theta)
	kappabf=sigmabf*graysaha
	kappabf=kappabf*(1.-np.exp(-h*c/(wav*1E-8*kerg*temp)))

	lwav=np.log10(wav)

	f0 = -2.2763 -1.6850*lwav +0.76661*lwav**2 -0.0533464*lwav**3
	f1 = 15.2827 -9.2846*lwav +1.99381*lwav**2 -0.142631*lwav**3
	f2 = ( -197.789 +190.266*lwav -67.9775*lwav**2 +10.6913*lwav**3 -0.625151*lwav**4 )

	ltheta=np.log10(theta)
	kappaff = 1E-26*elpress*10**(f0+f1*ltheta+f2*ltheta**2)
	
	return kappaff + kappabf

wav = np.linspace(3000,19000, 1000)
FALCexthmin = exthmin(wav,temp[np.where(h==0)][0],nel[np.where(h==0)][0])
print temp[np.where(h==0)][0]
print nel[np.where(h==0)][0]
print r'Maximum of exhtmin = %e at wavelength = %s angstrom$'%(np.max(FALCexthmin),
wav[np.where(FALCexthmin == np.max(FALCexthmin))][0]) 

plt.figure(25)
plt.title(r'H-minus bound-free + free-free extinction per H atom')
plt.plot(wav*1e-4,FALCexthmin*1e26/(nel[np.where(h==0)][0]*kerg*temp[np.where(h==0)][0]))	#Note the 10^-24 i units
plt.ylabel(r'$\kappa /P_e [10^{-26}$cm$^2$/dyn cm$^{-2}$]',size=14)
plt.xlabel(r'Wavelength [$\mu$m]')
plt.xlim([0,2])
plt.ylim([0,5])
#plt.xscale('log')
#plt.yscale('log')
#plt.show()

#Plot variation of H^- exctinction with height for wav = 0.5 \mu m
exthmin05 = exthmin(5000,temp,nel)
plt.figure(26)
plt.plot(h,exthmin05*(nhyd - nprot)*1e26)
plt.title(r'H-minus bf+ff extinction at $\lambda = 0.5 \mu m$')
plt.ylabel(r'log $(\kappa)$ [cm$^{-1}$]',size=14)
plt.xlabel(r'Height [km]')
plt.yscale('log')
plt.show()
