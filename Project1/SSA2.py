# importing useful libraries (you may need more)
import numpy as np # numerical package
import matplotlib.pyplot as plt # plotting package
from matplotlib import rc 
rc('font',**{'family':'serif'}) # This is for Latex writing

# Definition of constants and variables
k = 8.61734e-5 # Boltzmann constant in eV/deg
KeV = 8.61734e-5
kerg = 1.380658e-16	#Boltzmann constant in erg K 
kev = 8.6173e-5	#Boltzmann constant in eV/deg
h = 6.62607e-27	#Planck constant in erg s
elmass = 9.109390e-28	#electron mass in grams

# definition of functions
def partfunc_E(temp):
	"Parition function U_r"
	chiion = np.array([7,16,31,51])
	u = np.zeros(4)
	for r in range(4):
		for s in range(chiion[r]):
			u[r] = u[r] + np.exp(-s/(k*temp))
	return u 		#returns all the values of u array

def boltz_E(temp,r,s):
	"Boltzmann distribution n_r,s/N_r"
	u = partfunc_E(temp)
	relnrs = 1. / u[r-1]*np.exp(-(s-1)/(KeV*temp))
	return relnrs

def saha_E(temp,elpress,ionstage):
	keVT = kev*temp
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

#Main part calling the functions to do things


#Printing the population vs temperature for ground stage
temp = np.arange(0,30001,1000)
#print temp
pop = np.zeros((5,31))
for T in np.arange(1,31):	
	for r in np.arange(1,5):
		pop[r,T] = sahabolt_E(temp[T],131.,r,1)
labellst= ['ground stage','first ion stage','seconds ion stage','third ion stage']
print pop
plt.figure(0)
#plt.hold('on')
#ground state plot
for i in range(1,5):
	plt.plot(temp,pop[i,:],label=labellst[i-1])

plt.xlabel('temperature',size=14)
plt.ylabel('population',size=14)
plt.yscale('log')
plt.ylim([1e-3,1.1])
plt.legend(loc='best')
plt.show()


#test block
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

#compute partition function for T=5000,10000,20000
for temp in [5000,10000,20000]:
	print partfunc_E(temp)

for r in range(1,6):
	print saha_E(20000,1e3,r)
print "----"
for r in range(1,6):
	print saha_E(10000,1e3,r)
#Compute 
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
