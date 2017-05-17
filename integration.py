import numpy as np
import interpolation as inter

#Function which give the derivative
def interpolation_derivative_result(x, x_array, y_array, n, y_sec, i):
    if (i == n - 1):
        return y_array[n - 1]
    k = x_array[i + 1] - x_array[i]
    z = inter.zeta(x, x_array, i)

    inter_dx = (y_array[i + 1] - y_array[i])/k
    inter_dx = inter_dx + (k * y_sec[i] / 6) * (3 * (1 - z)**2 + 1)
    inter_dx = inter_dx + (k * y_sec[i+1] / 6) * (3 * z**2 - 1)
    return inter_dx

#Function which give the derivative of the interpolation 
def interpolation_derivative(x_array, y_array):
    n = x_array.size
    A = inter.create_matrix_A(x_array, n)
    B = inter.create_matrix_B(x_array, y_array, n)
    y_sec = inter.solve_linear_equations(A, B, n - 2)
    f = (lambda x: np.sqrt(1+ interpolation_derivative_result(x, x_array, y_array, n, y_sec, inter.dichotomic_search(x, x_array, 0, n))**2))
    return f

#------------------------------------------------------------------
#Functions implemented but not used because they are not efficient enough

def rectangle_method(f, n, a, b):
    h = (b - a) / n
    z = 0
    for i in range (n):
        z = z + f(a+i*h)
    return h*z


def pt_middle_method(f, n, a, b):
    h = (b - a) / n
    z = 0
    for i in range (n):
        z = z + f(((a+i*h) + (a+(i+1)*h)) /2) 
    return h*z

#------------------------------------------------------------------
#Functions implemented and used for the test in the file test_integration.py
#which give the integral using the different methods and take :
# -the function f to integrate
# -an interval [a, b]
# -the number of subdivisions

def trapezoidal_rule(f, n, a, b):
    h = (b - a) / n
    z = 0.5 * (f(a) + f(b))
    for i in range(1,n) :
        z = z + f(a+i*h)
    return h*z


def Simpson_rule(f, n, a, b) :
    h=(b - a) / n
    z=f(a) + f(b)
    for i in range(1, n) :
        z = z + 2 * f(a+i*h)
    for i in range(n) :
        z = z + 4 * f(a+(2*i+1)*h/2)
        
    return (h/6)*z


def Romberg_method(f, a, b, n):
    I = np.zeros((n, n))

    for k in range(0, n):    
        I[k, 0] = trapezoidal_rule(f, 2**k, a, b)
        
        for j in range(0, k):
            
            if (np.abs(I[k, j] - I[k-1, j-1]) < 10e-5):
                return I[k, j]
            
            I[k, j+1] = (4**(j+1) * I[k, j] - I[k-1, j]) / (4**(j+1) - 1)
            
        
    return I[n-1, n-1]


def Legendre(r, x):
    p0 = 1
    p1 = x
    for i in range (2, r+1):
        p = ((2*i-1)*p1 * x - (i-1) * p0)/ i
        p0 = p1
        p1 = p
    return p1

def Legendre_der(r, x):
    pr = Legendre(r, x)
    pr1 = Legendre(r + 1, x)
    return ((r + 1) / (1 - x**2)) * (x * pr - pr1)


#Returns the lambda-expression for the Legendre polynomial
def lambda_Legendre(r):
    return (lambda x: Legendre(r, x))

#Returns the lambda-expression for the derivate of Legendre polynomial
def lambda_Legendre_der(r):
    return (lambda x: Legendre_der(r, x))

#Returns the factorial of x
def fact(x):
    result = 1
    while (x>1) :
        result=result*x
        x=x-1
    return result


#Returns the coefficient of the i nth Legendre polynomial
def coefficient_Legendre(n):
    fin = n
    c = 1 / 2**n
    if (n%2 != 0):
        fin = n
        coef = np.zeros(n+1)
    else :
        coef = np.zeros(n+1)
    
    for i in range (fin // 2 + 1):
    
        f1 = fact(2*n - 2*i)
        f2 = fact(i)
        f3 = fact(n-i)
        f4 = fact(n - 2*i)
        f = ((-1)**i * f1) / (f2 * f3 * f4)
        coef[2*i] = c * f
        
    return (coef)

    
#Returns the weight list of each terms of a polynomial 
def weight_legendre(r, n):
    lld = lambda_Legendre_der(n)
    w = np.zeros(np.size(r))
    for i in range (np.size(r)):
        w[i] = 2 / ((1 - r[i]**2) * lld(r[i])**2)
    return w

def Gauss_Legendre(f, n, a, b):
    r =np.roots(coefficient_Legendre(n))
    w = weight_legendre(r, n)
    I = 0
    h = (b - a) / 2
    for i in range (n):
        I = I + w[i] * f(h * r[i] + ((a+b)/2))

    return h * I
