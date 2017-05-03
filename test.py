import numpy as np
import re #regexp
import matplotlib.pyplot as plt
import load_foil as ld
import interpolation as inter

def display(ex, ey, ix, iy):
    plt.plot(ex, ey, linewidth = 1.0)
    plt.plot(ix, iy, linewidth = 1.0)
    plt.show()

def display_interpole(ex, f_ex, ix, f_ix, n):
    ex_t = np.empty(0)
    ix_t = np.empty(0)
    for i in range(0, ex.size - 1):
        ex_tmp = np.arange(ex[i], ex[i + 1], (ex[i + 1] - ex[i])/n)
        for j in range(0, n):
            ex_t = np.append(ex_t, ex_tmp[j])

    for i in range(0, ix.size - 1):
        ix_tmp = np.arange(ix[i], ix[i + 1], (ix[i + 1] - ix[i])/n)
        for j in range(0, n):
            ix_t = np.append(ix_t, ix_tmp[j])

    ey_t = np.empty(ex_t.size)
    iy_t = np.empty(ix_t.size)
    for i in range(0, ex_t.size):
        ey_t[i] = f_ex(ex_t[i])
        
    for i in range(0, ix_t.size):
        iy_t[i] = f_ix(ix_t[i])

    plt.plot(ex_t, ey_t)
    plt.plot(ix_t, iy_t)
    plt.show()

(dim, ex,ey,ix,iy) = ld.load_foil("airfoils/b29root.dat")
display(ex, ey, ix, iy)
f = inter.interpolation(ex, ey)
g = inter.interpolation(ix, iy)
display_interpole(ex, f, ix, g, 10)
