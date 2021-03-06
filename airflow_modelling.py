import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.image as img
import interpolation as sp
import integration as it
import curves
from load_foil import load_foil

accuracy_x = 0.001
accuracy_y = 0.001
air_density = 1.225

#take a function en return an array. The array contains f(x) for each value of x_array
def function_to_array(f, x_array):
        y_array = np.empty(x_array.size)
        compt = 0
        for i in x_array:	
                y_array[compt] = f(i)
                compt = compt + 1
        return y_array

#color the image between two functions (function_min and function_max)
def coloring_image_part(image, function_min, function_max, x_min, x_max, y_min, value):
	#transform a pressure value to a color value.
	i = x_min
	while i < x_max:
                #find the index corresponding to this part of the image (accuracy_x/2 is used to be a the center of the "case")
		index_x = int(np.floor((i-x_min+(accuracy_x/2))/accuracy_x))
		max = function_max(i)
		j = function_min(i)
		while j < max:
			index_y = int(np.floor((j-y_min)/accuracy_y))
                        #color the image
			image[index_y][index_x] = [value, value, value]
			j = j + accuracy_y
		i = i + accuracy_x
		
#calcule the pressure for all the functions given
def compute_pressures(pressures, functions, x_array, y_min, y_max, index_start):
	x_min = x_array[0]
	x_max = x_array[x_array.size - 1]

	functions_n = len(functions)
	for i in range(0, functions_n, 1):
		f = it.interpolation_derivative(x_array , function_to_array(functions[i], x_array))
		length = it.Gauss_Legendre(f, x_array.size - 1, x_min, x_max)

		pressure = (air_density / 2) * (length**2)
		pressures[i + index_start] = pressure

#find the min and max of the pressures
def find_min_max_pressure(pressures):
	min_pressure = pressures[0]
	max_pressure = pressures[0]
	for pressure in pressures:
		if pressure < min_pressure:
			min_pressure = pressure
		if pressure > max_pressure:
			max_pressure = pressure
	return (min_pressure, max_pressure)

#calculate the value of pressure and color the upper or lower part of the image.
def creating_pressure_map(image, functions, x_array, y_min, y_max, up, pressures, min_pressure, max_pressure, index_start):
	x_min = x_array[0]
	x_max = x_array[x_array.size - 1]

	y_zero = (lambda x: 0)
	coloring_image_part(image, y_zero, functions[0], x_min, x_max, y_min, min_pressure)

	functions_n = len(functions)                        
	for i in range(0, functions_n, 1): 
		pressure = int(np.floor((pressures[i + index_start] - min_pressure) / (max_pressure-min_pressure) * 255))
		if i == functions_n - 1:
                        #color the extremities of the image.
			if up:
				y_max_f = (lambda x: y_max)
				coloring_image_part(image, functions[i], y_max_f, x_min, x_max, y_min, min_pressure)
			else:
				y_min_f = (lambda x: y_min)
				coloring_image_part(image, y_min_f, functions[i], x_min, x_max, y_min, min_pressure)
		else:
			if up:
				function_min = functions[i]
				function_max = functions[i + 1]
			else:
				function_min = functions[i + 1]
				function_max = functions[i]

			coloring_image_part(image, function_min, function_max, x_min, x_max, y_min, pressure)


#create the two sides of the image.
def create_image(airflow_up, airflow_up_n, airflow_down, airflow_down_n, x_array, y_min, y_max): 
	x_n = int(np.floor((x_array[x_array.size - 1] - x_array[0]) / accuracy_x))
	y_n = int(np.ceil((y_max - y_min) / accuracy_y))
	
	#image with color (hence the 3)
	image = np.empty([y_n, x_n, 3], dtype=np.uint8)
	
	pressures = np.zeros(airflow_up_n + airflow_down_n)
	
	compute_pressures(pressures, airflow_up, x_array, y_min, y_max, 0)

	compute_pressures(pressures, airflow_down, x_array, y_min, y_max, airflow_up_n)


	(min_pressure, max_pressure) = find_min_max_pressure(pressures)
	
	#create the image.
	creating_pressure_map(image, airflow_up, x_array, y_min, y_max, 1, pressures, min_pressure, max_pressure, 0)
	creating_pressure_map(image, airflow_down, x_array, y_min, 0, 0, pressures, min_pressure, max_pressure, airflow_up_n)
	return image

#create the pressure image
def pressure_mapping(input):
	# Maps the pressure around a wing given by input (a .dat file)"
        (dim, ex, ey, ix, iy) = load_foil(input)

        #value the farest of the center.
        h_max = np.max(ey)
        h_min= np.min(iy)
	
        #interpolation function of the wing
        f_up = sp.interpolation(ex,ey)
        f_low = sp.interpolation(ix,iy)

	#number of air slices for each part of the wing.
        n_curves = 10

        #airflow above the wing
        f_up_lambdas = curves.create_curves(f_up, n_curves, h_max)

        #airflow below the wing
        f_low_lambdas = curves.create_curves(f_low, n_curves, h_min)

        #creation of the simulation of the airflow on the wing
        image = create_image(f_up_lambdas, len(f_up_lambdas), f_low_lambdas, len(f_low_lambdas), ex, 3*h_min, 3*h_max)
        lum_imge = image[:,:,0]
        my_cmap = plt.cm.get_cmap('hot')
        display = plt.imshow(lum_imge, cmap=my_cmap, vmin = 0, vmax = 255, interpolation='nearest' ,origin='lower')
        print("Pressure mapping done")
        plt.show()
        
pressure_mapping("airfoils/b29root.dat")
