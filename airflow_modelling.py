import numpy as np
import pylab as pl
import interpolation as sp
import integration as it
import curves
from load_foil import load_foil

def f_lambda(f, x, h):
	# Returns the function f_lambda, with h = hmin or hmax
	return lambda x: (1 - x) * f(x) + 3 * x * h

def show_colored(f, a, b, colour = (0,0,1), width = 4):
	# Create the graph of a function f on the interval [a,b], in color"
	v = []
	fv = []
	steps = 256.
	for i in np.arange(a, b + (b - a)/steps, (b - a)/steps):
		fv.append(f(i))
	pl.plot(np.arange(a, b + (b - a)/steps, (b - a)/steps), fv, color = colour, linewidth = width)

def get_color(Pmin, Pmax, P):
	# For a given pressure, gives the right color to use"
	color = ((1,1,1), (1,1,0), (1,204/255.,102/255.), (1,153/255.,0), (1,0,0), (204/255.,0,0), (153/255.,0,0), (102/255.,0,0), (51/255.,0,0), (0,0,0))
	pitch = (Pmax**2 - Pmin**2) / 100 # 100 = square of the number of colors
	P = P**2 - Pmin**2

	if (0<=P<=pitch):
		return color[9]
	elif (pitch < P <= 2* pitch):
		return color[8]
	elif (2 * pitch < P <= 3* pitch):
		return color[7]
	elif (3 * pitch < P <= 4* pitch):
		return color[6]
	elif (4 * pitch < P <= 5 * pitch):
		return color[5]
	elif (5* pitch < P <= 6 * pitch):
		return color[4]
	elif (6 * pitch < P <= 7 * pitch):
		return color[3]
	elif (7 * pitch < P <= 8 * pitch):
		return color[2]
	elif (8 * pitch < P <= 9 * pitch):
		return color[1]
	else:
		return color[0]

def pressure_mapping(input, pitch = 0.01):
	# Maps the pressure around a wing given by input (a .dat file)"
	pl.clf()
	(dim, ex, ey, ix, iy) = load_foil(input) # Outside function
	h_max = np.max(ey)
	h_min = np.min(iy)

	f_up = sp.interpolation(ex,ey) # Outside function
	f_low = sp.interpolation(ix,iy) # Outside function

	n_curves = 10
	
	f_up_lambdas = curves.create_curves(f_up, n_curves, h_max)
	f_low_lambdas = curves.create_curves(f_low, n_curves, h_min)
	
	precision_length = np.arange(0,1,0.01)

	Pmax = it.Gauss_Legendre(f_up, 10, ex[0], ex[ex.size - 1])
	Pmin = it.Gauss_Legendre(f_low, 10, ix[0], ix[ix.size - 1])
	#Pmax = it.length(sp.interpolation(ex,ey),0.,1.,100) # Two outside functions
	#Pmin = it.length(sp.interpolation(np.arange(0,1,0.01), map(f_lambda(fup,1,hmax),precision_length)), 0., 1., 100) # Three outside functions
	#for i in np.arange(1e-2 , 2 , pitch): # over the extrados, 1e-2 is to avoid superposition of the wing and the pressure
		#f = f_lambda(f_up , i  ,h_max)
	for f in f_up_lambdas: 
		v = it.Gauss_Legendre(f, 10, ex[0], ex[ex.size - 1])
		#v = it.Gauss_Legendre(sp.interpolation(np.arange(0,1,0.01), map(f,precision_length)) ,0.,1.,100)
		show_colored(f , ex[0],ex[ex.size - 1], get_color(Pmin, Pmax, v**2))
		#for i in np.arange(1e-2, 1 , pitch): # under the intrados
	#f = f_lambda(f_low , i ,h_min)
	for f in f_low _lambdas:
		v = it.Gauss_Legendre(f, 10, ix[0], ix[ix.size - 1])		
		#v = it.Gauss_Legendre(sp.interpolation(np.arange(0,1,0.01), map(f,np.arange(0,1,0.01))) ,0.,1.,100)
		show_colored(f, ix[0], ix[ix.size - 1], get_color(Pmin,Pmax, v**2))

	show_colored(f_up, ex[0], ex[len(ex) - 1],'g',2)
	show_colored(f_low, ix[0], ix[len(ix) -1],'g',2)

	print("Pressure mapping done")
	pl.show()

pressure_mapping("airfoils/b29root.dat")
