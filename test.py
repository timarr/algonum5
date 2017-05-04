import numpy as np
import re #regexp
import matplotlib.pyplot as plt
import load_foil as ld
import interpolation as inter
import curves

def display(ex, ey, ix, iy):
    plt.plot(ex, ey, linewidth = 1.0)
    plt.plot(ix, iy, linewidth = 1.0)
    plt.show()

def plot_funct(x, f_x, n):
    x_t = np.empty(0)
    for i in range(0, x.size - 1):
        x_tmp = np.arange(x[i], x[i + 1], (x[i + 1] - x[i])/n)
        for j in range(0, n):
            x_t = np.append(x_t, x_tmp[j])

    y = np.empty(x_t.size)
    for i in range(0, x_t.size):
        y[i] = f_x(x_t[i])

    plt.plot(x_t, y) 
    

def display_interpole(ex, f_ex, ix, f_ix, n):
    plot_funct(ex, f_ex, n)
    plot_funct(ix, f_ix, n)
    plt.show()



def display_slices(ex, f_lambdas, ix, g_lambdas, n):
    for f in f_lambdas:
        plot_funct(ex, f, n)
        
    for g in g_lambdas:
        plot_funct(ix, g, n)
    plt.show()


(dim, ex,ey,ix,iy) = ld.load_foil("airfoils/b29root.dat")
#display(ex, ey, ix, iy)

number_points = 10
f = inter.interpolation(ex, ey)
g = inter.interpolation(ix, iy)
#display_interpole(ex, f, ix, g, number_points)

n_curves = 10
h_max_e = curves.find_farrest_point(ey)
h_max_i = curves.find_farrest_point(iy)
f_lambdas = curves.create_curves(f, n_curves, h_max_e)
g_lambdas = curves.create_curves(g, n_curves, h_max_i)
display_slices(ex, f_lambdas, ix, g_lambdas, number_points)
