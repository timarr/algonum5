import numpy as np

def find_farrest_point(y_array):
    max_p = y_array[0]
    for i in y_array:
        if (i**2 > max_p**2):
            max_p = i
    return max_p

def create_curves(function, n_curves, h_max):
    curves = []
    lambdas = np.linspace(0., 1., num=n_curves)
    for i in range(0, lambdas.size):
        y = create_lambda(function, lambdas[i], h_max)
        curves.append(y)
    return curves

def create_lambda(function, value, h_max):
    return (lambda x: (1 - value) * function(x) + value * 3 * h_max)
