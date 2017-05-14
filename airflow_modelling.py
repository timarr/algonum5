import numpy as np
import scipy.integrate as sc
import scipy.misc as sm
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.image as img
import interpolation as sp
import integration as it
import curves
from load_foil import load_foil

accuracy_x = 0.001
accuracy_y = 0.001
static_pressure = 1
air_density = 1.225

def function_to_array(f, x_array):
        y_array = np.empty(x_array.size)
        for i in x_array:
                np.append(y_array, f(i))
        return y_array

def function_to_derivative(f, x_array):
        h = (x_array[x_array.size - 1] - x_array[0])/x_array.size
        y_array = np.empty(x_array.size)
        print("derivative")
        for i in x_array:
                y = sm.derivative(f, i, dx=h)
                np.append(y_array, y)
        return y_array
                

#color the image between two functions (function_min and function_max)
def coloring_image_part(image, function_min, function_max, x_min, x_max, y_min, value):
        i = x_min
        while i < x_max:
                #find the index corresponding to this part of the image (accuracy_x/2 is used to be a the center of the "case")
                index_x = np.floor((i-x_min+(accuracy_x/2))/accuracy_x)
                max = function_max(i)
                j = function_min(i)
                while j < max:
                        index_y = np.floor((j-y_min)/accuracy_y)
                        print(value)
                        value_i = np.floor(value)
                        #color the image
                        image[index_y][index_x] = [value_i, value_i, value_i]
                        j = j + accuracy_y
                i = i + accuracy_x

#calculate the value of pressure and color the upper or lower part of the image.
def creating_pressure_map(image, functions, x_array, y_min, y_max, poly_n, up):
        x_min = x_array[0]
        x_max = x_array[x_array.size - 1]
        #inside the wing
        y_zero = (lambda x: 0)
        coloring_image_part(image, y_zero, functions[0], x_min, x_max, y_min, 0)

        #number of functions
        functions_n = len(functions)
        for i in range(0, functions_n, 1):
                #length of the wing
                #f = it.interpolation_derivative(x_array , function_to_array(functions[i], x_array))
                #length = it.rectangle_method(f, poly_n, x_min, x_max)
                f = function_to_derivative(functions[i], x_array)
                f_prim = sp.interpolation(x_array, f)
                g = (lambda x: np.sqrt(1 + pow(f_prim(x),2)))
                length = it.Gauss_Legendre(g, x_array.size, x_min, x_max)
                dynamic_pressure = (air_density / 2) * (length**2)
                print(i)
                print(length)
                print(dynamic_pressure)
                print("--------")
                #pressure apply by the airflow
                pressure = dynamic_pressure
                if i == functions_n - 1:
                        #color the extremities of the image.
                        if up:
                                y_max_f = (lambda x: y_max)
                                coloring_image_part(image, functions[i], y_max_f, x_min, x_max, y_min, pressure)
                        else:
                                y_min_f = (lambda x: y_min)
                                coloring_image_part(image, y_min_f, functions[i], x_min, x_max, y_min, pressure)
                else:
                        if up:
                                function_min = functions[i]
                                function_max = functions[i + 1]
                        else:
                                function_min = functions[i + 1]
                                function_max = functions[i]

                        coloring_image_part(image, function_min, function_max, x_min, x_max, y_min, length)

#create the two sides of the image.
def create_image(airflow_up, airflow_up_n, airflow_down, airflow_down_n, x_array, y_min, y_max): 
        #size of the image
        x_n = np.floor((x_array[x_array.size - 1] - x_array[0]) / accuracy_x)
        y_n = np.ceil((y_max - y_min) / accuracy_y)
        #image of RGB (3)
        image = np.zeros([y_n, x_n, 3], dtype=np.uint8)
        creating_pressure_map(image, airflow_up, x_array, y_min, y_max, airflow_up_n, 1)
        creating_pressure_map(image, airflow_down, x_array, y_min, 0, airflow_down_n, 0)
        return image

def pressure_mapping(input, pitch = 0.01):
	# Maps the pressure around a wing given by input (a .dat file)"
        (dim, ex, ey, ix, iy) = load_foil(input)
        #value the farest of the center.
        h_max = np.max(ey)
        h_min= np.min(iy)
        #interpolation function of the wing
        f_up = sp.interpolation(ex,ey)
        
        f_low = sp.interpolation(ix,iy)

        n_curves = 10

        #airflow above the wing
        f_up_lambdas = curves.create_curves(f_up, n_curves, h_max)

        #airflow below the wing
        f_low_lambdas = curves.create_curves(f_low, n_curves, h_min)

        #creation of the simulation of the airflow on the wing
        image = create_image(f_up_lambdas, 10, f_low_lambdas, 10, ex, 3*h_min, 3*h_max)
        lum_imge = image[:,:,0]
        my_cmap = plt.cm.get_cmap('hot')
        display = plt.imshow(lum_imge, cmap=my_cmap, vmin = 0, vmax = 255, interpolation='nearest' ,origin='lower')
        print("Pressure mapping done")
        plt.show()
        
pressure_mapping("airfoils/b29root.dat")
