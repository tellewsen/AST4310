# importing useful libraries (you may need more)
import numpy as np # numerical package
import matplotlib.pyplot as plt # plotting package
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

# Definition of constants and variables
keV	= 8.61734e-5 # Boltzmann constant in eV/deg
kerg	= 1.380658e-16	#Boltzmann constant in erg K 
h	= 6.62607e-27	#Planck constant in erg s
c	= 2.99792e10	#Speed of light in cm/s
elmass	= 9.109390e-28	#electron mass in grams

# definition of functions
def partfunc_E(temp):
	"Parition function U_r"
	keVT = keV*temp
	chiion = np.array([7,16,31,51])
	u = np.zeros(4)
	for r in range(4):
		for s in range(chiion[r]):
			u[r] = u[r] + np.exp(-s/(keVT))
	return u 		#returns all the values of u array

def boltz_E(temp,r,s):
	"Boltzmann distribution n_r,s/N_r"
	u = partfunc_E(temp)
	relnrs = 1. / u[r-1]*np.exp(-(s-1)/(keV*temp))
	return relnrs

def saha_E(temp,elpress,ionstage):
	keVT = keV*temp
	kergT =kerg*temp
	eldens = elpress/kergT
	chiion = np.array([7,16,31,51])
	u = partfunc_E(temp)
	u = np.append(u,2) #append element to array
	sahaconst = (2* np.pi *elmass*kergT/(h**2))**(3./2)*2/eldens
	nstage = np.zeros(5)
	nstage[0] = 1.
	for r in range(4):
		nstage[r+1] =nstage[r]*sahaconst*u[r+1]/u[r]*np.exp(-chiion[r]/keVT) 	
	ntotal = np.sum(nstage)
	nstagerel = nstage/ntotal
	return nstagerel[ionstage-1]

def sahabolt_E(temp, elpress, ion, level):
	return saha_E(temp, elpress, ion)*boltz_E(temp, ion, level)

def sahabolt_H(temp,elpress,level):
	keVT = keV*temp
	kergT =kerg*temp
	eldens = elpress/kergT
	#energy levels and weights for hydrogen
	nrlevels = 100		#Pick a reasonable cutoff value for partition function
	g = np.zeros((2,nrlevels))	#declaration weights
	chiexc = np.zeros((2,nrlevels))	#declaration excitation energies

	for s in range(nrlevels):
		g[0,s] = 2.*(s+1.)**2			#statistical weights
		chiexc[0,s] = 13.598*(1.-1./(s+1.)**2.)	#excitation weights
	g[1,0] = 1.
	chiexc[1,0] = 0.

	#Partition functions
	u = np.zeros([2])
	for s in range(nrlevels):
		u[0] = u[0] + g[0,s]*np.exp(-chiexc[0,s]/keVT)
	u[1] = g[1,0]

	#Saha
	sahaconst = (2.* np.pi *elmass*kergT/(h*h))**(1.5)*2/eldens
	nstage = np.zeros(2)
	nstage[0] = 1.
	nstage[1] = nstage[0] * sahaconst*u[1]/u[0] *np.exp(-13.598/keVT)
	ntotal = np.sum(nstage)	#sum both stages = total hydrogen density
	
	#Boltzmann
	nlevel = nstage[0]*g[0,level-1]/u[0]*np.exp(-chiexc[0,level-1]/keVT)
	nlevelrel = nlevel/ntotal	#Fraction of total hydrogen density

	
	#Test block
	"""
	for s in range(6):
		print s+1, g[0,s],chiexc[0,s],g[0,s]*np.exp(-chiexc[0,s]/keVT)
	for s in range(0,nrlevels,10):
		print s+1, g[0,s],chiexc[0,s],g[0,s]*np.exp(-chiexc[0,s]/keVT)
	"""	


	return nlevelrel

def partfunc_Ca(temp):
	"Parition function U_r"
	keVT = keV*temp
	chiion = np.array([6.113, 11.871, 50.91, 67.15])
	u = np.zeros(4)
	for r in range(4):
		for s in range(int(chiion[r])):
			u[r] = u[r] + np.exp(-s/(keVT))
	return u 		#returns all the values of u array

def boltz_Ca(temp,r,s):
	"Boltzmann distribution n_r,s/N_r"
	keVT = keV*temp
	u = partfunc_Ca(temp)
	relnrs = 1. / u[r-1]*np.exp(-(s-1)/(keVT))
	return relnrs

def saha_Ca(temp,elpress,ionstage):
	keVT = keV*temp
	kergT =kerg*temp
	eldens = elpress/kergT
	chiion = np.array([6.113, 11.871, 50.91, 67.15])
	u = partfunc_Ca(temp)
	u = np.append(u,2) #append element to array
	sahaconst = (2.*np.pi *elmass*kergT/(h*h))**(3./2)*2./eldens
	nstage = np.zeros(5)
	nstage[0] = 1.
	for r in range(4):
		nstage[r+1] =nstage[r]*sahaconst*u[r+1]/u[r]*np.exp(-chiion[r]/keVT) 	
	ntotal = np.sum(nstage)
	nstagerel = nstage/ntotal
	return nstagerel[ionstage-1]

def sahabolt_Ca(temp, elpress, ionstage, level):
	return saha_Ca(temp, elpress, ionstage)*boltz_Ca(temp, ionstage, level)

#Main part calling the functions to do things

#Part2

#Compute partition function for T=5000,10000,20000
"""
for temp in [5000,10000,20000]:
	print partfunc_E(temp)
"""
"""
#Compute distribution for T=5000,
distribution_5k = np.zeros(11)
distribution_10k = np.zeros(11)
distribution_20k = np.zeros(11)
for s in range(1,11):
	distribution_5k[s-1] = boltz_E(5000,1,s)
	distribution_10k[s-1] = boltz_E(10000,1,s)
	distribution_20k[s-1] = boltz_E(20000,1,s)
print distribution_5k
print distribution_10k
print distribution_20k
"""
#Compute saha_E for different temp
"""
for r in range(1,6):
	print saha_E(20000,1e3,r)
print "----"
for r in range(1,6):
	print saha_E(10000,1e3,r)
"""

#Compute Saha_bolt_E for different temperatures
"""
for s in range(1,6):
	print sahabolt_E(5000, 1e3, 1, s)
print "----"
for s in range(1,6):
	print sahabolt_E(20000, 1e3, 1, s)
print "----"
for s in range(1,6):
	print sahabolt_E(10000, 1e3, 2, s)
print "----"
for s in range(1,6):
	print sahabolt_E(20000, 1e3, 4, s)
"""
#Plotting the population vs temperature for s=1 and P_e = 131
"""
temp = np.arange(0,30001,1000)
#print temp
pop = np.zeros((5,31))
for T in np.arange(1,31):	
	for r in np.arange(1,5):
		pop[r,T] = sahabolt_E(temp[T],131.,r,1)
labellst= ['ground stage','first ion stage','seconds ion stage','third ion stage']
#print pop
#ground state plot
plt.figure(0)
for i in range(1,5):
	plt.plot(temp,pop[i,:],label=labellst[i-1])

plt.xlabel('temperature',size=14)
plt.ylabel('population',size=14)
plt.yscale('log')
plt.ylim([1e-3,1.1])
plt.legend(loc='best')
plt.title('Population vs Temperature for s=1')
plt.show()
"""
#Plotting the population vs temperature for s=2
"""
temp = np.arange(0,30001,1000)
#print temp
pop = np.zeros((5,31))
for T in np.arange(1,31):	
	for r in np.arange(1,5):
		pop[r,T] = sahabolt_E(temp[T],131.,r,2)
labellst= ['ground stage','first ion stage','seconds ion stage','third ion stage']
#print pop
#ground state plot
plt.figure(1)
for i in range(1,5):
	plt.plot(temp,pop[i,:],label=labellst[i-1])

plt.xlabel('temperature',size=14)
plt.ylabel('population',size=14)
plt.yscale('log')
plt.ylim([1e-3,1.1])
plt.legend(loc='best')
plt.title('Population vs Temperature for s=2')
plt.show()
"""
#Plotting the population vs temperature for s=4
"""
temp = np.arange(0,30001,1000)
#print temp
pop = np.zeros((5,31))
for T in np.arange(1,31):	
	for r in np.arange(1,5):
		pop[r,T] = sahabolt_E(temp[T],131.,r,4)
labellst= ['ground stage','first ion stage','seconds ion stage','third ion stage']
#print pop
#ground state plot
plt.figure(2)
for i in range(1,5):
	plt.plot(temp,pop[i,:],label=labellst[i-1])

plt.xlabel('temperature',size=14)
plt.ylabel('population',size=14)
plt.yscale('log')
plt.ylim([1e-3,1.1])
plt.legend(loc='best')
plt.title('Population vs Temperature for s=4')
plt.show()
"""
#Printing hydrogen levels
"""
sahabolt_H(5000,1e2,1) this works as it should
"""

#Solar Ca+K versus Ha:line strength
"""
temp = np.arange(1000,20001,100)
CaH = np.zeros(temp.shape)
Caabund = 2e-6
for i in range(0,len(temp)):
	NCa = sahabolt_Ca(temp[i],1e2,2,1)
	NH  = sahabolt_H(temp[i],1e2,2)
	CaH[i] = NCa*Caabund/NH

plt.plot(temp,CaH,label=r'strength ratio Ca$^+$K /H$\alpha$')
plt.yscale('log')
plt.xlabel(r'temperature $T / K$',size=14)
plt.ylabel(r'Ca II K / H$\alpha$',size=14)
plt.legend(fontsize=14)
plt.title('Ca/H Ratio versus temperature')
plt.show()
print 'Ca/H ratio at 5000 K = ', CaH[np.argwhere(temp==5000)][0][0]

"""
#solar Ca+K versus Ha: temperature sensitivity
"""
temp = np.arange(2000,12001,100)
dNCadT = np.zeros(temp.shape)
dNHdT = np.zeros(temp.shape)
dT = 1.
for i in range(101):
	NCa = sahabolt_Ca(temp[i],1e2,2,1)
	NCa2 = sahabolt_Ca(temp[i]-dT,1e2,2,1)
	dNCadT[i] = (NCa - NCa2)/(dT*NCa)
	NH = sahabolt_H(temp[i],1e2,2)
	NH2 = sahabolt_H(temp[i]-dT,1e2,2)
	dNHdT[i] = (NH-NH2)/(dT*NH)
plt.figure(0)
plt.plot(temp,np.abs(dNHdT),label=r'H')
plt.plot(temp,np.abs(dNCadT),ls='--',label=r'Ca$^+$K')
plt.yscale('log')

plt.title('Relative population changes')
plt.ylim(1e-9,1e1)

NCa = np.zeros(temp.shape)
NH = np.zeros(temp.shape)
for i in range (101):
	NCa[i] = sahabolt_Ca(temp[i],1e2,2,1)
	NH[i] = sahabolt_H(temp[i],1e2,2)
plt.plot(temp,NH/np.amax(NH), label = 'rel. pop H')
plt.plot(temp,NCa/np.amax(NCa),ls='--',label = r'rel. pop Ca$^+$K')

plt.xlabel(r'temperature $T/K$',size=14)
plt.ylabel(r'$\left| \left( \Delta n(r,s) / \Delta T\right) /  n(r,s) \right|$',size=20)
plt.legend(loc=4,fontsize=12)	
plt.show()
"""

#Hot stars vs cold stars
"""
"Find at which temperature the hydrogen in stellar photospheres with P_e = 100 is about 50% ionized."
for T in np.arange(7000,10001,100):
	print T,sahabolt_H(T,1e2,1)
temp = np.arange(6000,13001,1e2)
nH = np.zeros(temp.shape)
for i in range(len(temp)):
	nH[i] = sahabolt_H(temp[i],1e2,1)

plt.plot(temp,nH)
plt.xlabel(r'temperature $T/K$',size=14)
plt.ylabel(r'neutral hydrogen fraction',size=14)
plt.title('Fraction of neutral hydrogen in stellar photospheres')
plt.show()
"""
