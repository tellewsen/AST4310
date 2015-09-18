# importing useful libraries
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

#Part3
def planck(temp,wav):
	B = 2.*h*c*c/(wav**5) / (np.exp( (h*c) / (wav*kerg*temp) ) - 1 )
	return B

#Planck function test ok
#print planck(5000,5000e-8)

wav = np.arange(1000,20801,200)
b = np.zeros(wav.shape)

plt.xlabel(r'wavelength $\lambda / \AA$',size=14)
plt.ylabel(r'Planck function',size=14)
plt.title('Radiation intensity emitted')
for T in np.arange(5000,8001,200):
	b[:] = planck(T,wav[:]*1e-8)
	plt.plot(wav,b,'-')
plt.yscale('log')
#plt.xscale('log')
plt.xlim(0,20800)
plt.show()
