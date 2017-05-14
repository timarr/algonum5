import numpy as np
from scipy import interpolate
 
# x : the points
# y : the y of each point
# n : the amount of points
# yp1 : f(x1)
# ypn : f(xn)
def spline(x, y, n, yp1, ypn):
    u = np.zeros((1, n))
    i = 2
    y2 = np.zeros((1, n)) #Will contain the derivatives of the points
    if (yp1 > 0.99e30):
        u[0][0] = 0
        y2[0][0] = 0.0
    else:
        y2[0][0] = -0.5
        u[0][0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0]))-yp1

    for i in range (1, n-1) :
        sig = (x[i]-x[i-1])/(x[i+1]-x[i-1])
        p = sig*y2[0][i-1]+2
        y2[0][i] = (sig-1)/p
        u[0][i] = ((y[i+1] - y[i]) / (x[i+1] - x[i])) - ((y[i] - y[i-1]) / (x[i] - x[i-1]))
        u[0][i] = (6.0*u[0][i]/(x[i+1] - x[i-1])-sig*u[0][i-1])/p
    if (ypn > 0.99e30):
        qn, u = 0.0, 0.0
    else :
        qn = 0.5
        un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]))
    y2[0][n-1] = (un - qn*u[0][n-2])/(qn*y2[0][n-2]+1)

    for k in range (n-2, 0, -1):
        y2[0][k] = y2[0][k]*y2[0][k+1]+u[0][k]
    return y2

#Computes the array xa[1..n] and ya[1..n]
def splint(xa, ya, y2a, n, x):
    klo = 0
    khi = n-1
    while (khi-klo > 1):
        k = (khi + klo) / 2
        if (xa[k] > x):
            khi = k
        else :
            klo = k
    h=xa[khi] - xa[klo]
    if (h == 0):
        print("error")
    a= (xa[khi]-x)/h
    b = (x - xa[klo])/h
    y = a*ya[klo] + b*ya[khi] + (((a*a*a-a)*y2a[0][klo])+(b*b*b-b)*y2a[0][khi])*(h*h)/6.0
    return y


