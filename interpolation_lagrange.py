import numpy as np
from scipy import interpolate
 
def interpolate_lagrange(dx,dy,data,x):
	def b(j,xi):
		v = 1.0
		for k in xrange(len(data)):
			if k != j:
				v *= (xi-dx[k]) / (dx[j]-dx[k])
		return v
	for i,xi in enumerate(dx):
		for j in xrange(len(data)):
			x += dy[j]*b(j,xi)
	return x
	
def lagrange_step(data, x):
	
	#unzipping
    dx = [d[0] for d in data] #stores into dx the first of the each couple
    dy = [d[1] for d in data] #same for the second of each couple
    return (lambda x : interpolate_lagrange(dx,dy,data,x))


