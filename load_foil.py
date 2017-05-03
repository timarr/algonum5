import numpy as np;
import re; # regexp
import matplotlib.pyplot as ma;

################################################################
# Airfoil : load profile of a wing
#
# Reads a file whose lines contain coordinates of points,
# separated by an empty line.
# Every line not containing a couple of floats is discarded. 
# Returns a couple constitued of the list of points of the
# extrados and the intrados. 
def load_foil(file):
    f = open(file, 'r')
    matchline = lambda line: re.match(r"\s*([\d\.-]+)\s*([\d\.-]+)", line)
    extra  = [];    intra = []
    rextra = False; rintra = False; rheader = False;
    for line in f:
        m = matchline(line)
        if (m == None):
            if not(rheader):
                rheader = True
            elif not(rextra):
                rextra = True
            elif not(rintra):
                rintra = True
            continue;
        if (m != None) and rheader and not(rextra) and not(rintra):
            dim = np.array(map(lambda t: float(t), m.groups()))
            continue
        if (m != None) and rheader and rextra and not(rintra):
            extra.append(m.groups())
            continue;
        if (m != None) and rheader and rextra and rintra:
            intra.append(m.groups())
            continue;
    i = 0
    ex = np.empty(len(extra))
    ey = np.empty(len(extra))
    for t in extra:
        ex[i] = t[0]
        ey[i] = t[1]
        i = i + 1

    i = 0
    ix = np.empty(len(intra))
    iy = np.empty(len(intra))
    for t in intra:
        ix[i] = t[0]
        iy[i] = t[1]
        i = i + 1
    return(dim,ex,ey,ix,iy)

