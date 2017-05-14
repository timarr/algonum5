import numpy as np
import re #regexp
import matplotlib.pyplot as plt
import load_foil as ld
import interpolation as inter
import curves
import interpolation_lagrange as lagrange

#Displays the plot describing the original wing
def display(ex, ey, ix, iy):
    plt.plot(ex, ey, linewidth = 1.0)
    plt.plot(ix, iy, linewidth = 1.0)
    plt.title("Courbes originales decrivant la forme de l'aile du b-29")
    plt.show()

#Generates the plot x->f(x)
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

#Displays the interpolated points of given curves    
def display_interpole(ex, f_ex, ix, f_ix, n):
    plot_funct(ex, f_ex, n)
    plot_funct(ix, f_ix, n)
    plt.title("Courbes interpolees (tri-diagonale matrice methode)")
    plt.show()
    
	
#Displays all the curves disposed around the wing modeling the airflow
def display_slices(ex, f_lambdas, ix, g_lambdas, n):
    for f in f_lambdas:
        plot_funct(ex, f, n)
        
    for g in g_lambdas:
        plot_funct(ix, g, n)
    plt.title("Courbes représentant la modification d'un courant d'air par l'aile")
    plt.show()

try:
    number_points=int(input('Entrez un nombre de points avec lequel vous voulez interpoler :'))
except ValueError:
    print ("Ceci n'est pas un nombre !")

    
print("Courbes originales decrivant la forme de l'aile du b-29")
(dim, ex,ey,ix,iy) = ld.load_foil("airfoils/b29root.dat")
display(ex, ey, ix, iy)

    
f = inter.interpolation(ex, ey)
g = inter.interpolation(ix, iy)
display_interpole(ex, f, ix, g, number_points)


x = np.linspace(0, 1, number_points)
y = []
splinval = lagrange.spline(ex, ey, len(ex), ey[0], ey[-1])
for i in range(len(x)):
	y.append(lagrange.splint(ex, ey, splinval, len(ex), x[i]))
plt.plot(x, y,'b')
y = []
splinval = lagrange.spline(ix, iy, len(ix), iy[0], iy[-1])
for i in range(len(x)):
	y.append(lagrange.splint(ix, iy, splinval, len(ex), x[i]))
plt.plot(x, y,'r')
plt.title("Courbes interpolees (methode de Lagrange)")
plt.show()




n_curves = 10
h_max_e = curves.find_farrest_point(ey)
h_max_i = curves.find_farrest_point(iy)
f_lambdas = curves.create_curves(f, n_curves, h_max_e)
g_lambdas = curves.create_curves(g, n_curves, h_max_i)
display_slices(ex, f_lambdas, ix, g_lambdas, number_points)
