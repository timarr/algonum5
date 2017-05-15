import numpy as np
import matplotlib.pyplot as plt
from integration import *

#Returns the integral of sqrt(x) between a and b
def integrale_racine_x(a, b):
    return ((2/3)*b**(3/2) - (2/3)*a**(3/2))

#Test the integration's methods whith the function sqrt(x)
def test_convergence(a, b, n):
    racine_x = lambda x: np.sqrt(x)
    
    t = np.zeros(n)
    s = np.zeros(n)
    r = np.zeros(n)
    g = np.zeros(n)

    for i in range (1 ,n+1):

        t[i-1] = np.linalg.norm(trapezoidal_rule(racine_x, i, a, b) - (integrale_racine_x(a, b))) / np.linalg.norm(integrale_racine_x(a, b))
        
        s[i-1] = np.linalg.norm(Simpson_rule(racine_x, i, a, b) - (integrale_racine_x(a, b))) / np.linalg.norm(integrale_racine_x(a, b))

        r[i-1] = np.linalg.norm(Romberg_method(racine_x, a, b, i) - (integrale_racine_x(a, b))) / np.linalg.norm(integrale_racine_x(a, b))

        g[i-1] = np.linalg.norm(Gauss_Legendre(racine_x, i, a, b) - (integrale_racine_x(a, b))) / np.linalg.norm(integrale_racine_x(a, b))
        
    plt.plot(t, label='Trapezoidal Rule')
    plt.plot(s, label='Simpson Rule')
    plt.plot(r, label='Romberg Method')
    plt.plot(g, label='Gauss-Legendre')

    plt.title("Convergence of the different methods of integrals - Function square root  ")
    plt.legend()
    plt.show()


test_convergence(0, 10, 10)   
