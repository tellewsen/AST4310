from numpy import *

"""
sums 1D array ARR (but IDL's total is faster and more general)
"""
arr = array([2])

def addup(arr):
	sumarr=0
	for i in range(len(arr)):
	    sumarr +=arr[i]
	return sumarr

print addup(arr)

