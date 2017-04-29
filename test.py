import numpy as np
import re #regexp
import matplotlib.pyplot as plt
import load_foil as ld


def display(ex, ey, ix, iy):
    plt.plot(ex, ey, linewidth = 1.0)
    plt.plot(ix, iy, linewidth = 1.0)
    plt.show()


(dim, ex,ey,ix,iy) = ld.load_foil("airfoils/b29root.dat")
display(ex, ey, ix, iy)

