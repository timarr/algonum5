import numpy as np

def create_matrix_A (x_array, n):
    n = x_array.size
    A = np.zeros((n - 2, n - 2))
    k = np.empty(n - 1)
    for i in range (0, n - 1):
        k[i] = x_array[i + 1] - x_array[i]
        
    for i in range (1, n - 2) :
        A[i - 1][i - 1] = (k[i] + k[i - 1]) / 3
        A[i][i - 1] = k[i - 1] / 6
        A[i - 1][i] = k[i] / 6 
    A[n - 3][n - 3] = k[n - 3] + k[n - 2]
    return A
        
def create_matrix_B (x_array, y_array, n):
    B = np.empty(n - 2) 

    k = np.empty(n - 1)
    for i in range (0, n - 1):
        k[i] = x_array[i + 1] - x_array[i]
        
    for i in range (1, n - 1):
        B[i - 1] = ((y_array[i + 1] - y_array[i]) / k[i]) - ((y_array[i] - y_array[i - 1]) / k[i - 1])

    return B


def solve_linear_equations(A, B, n):
    y_sec = np.empty(n + 2)
    y_sec[0] = 0
    y_sec[n + 1] = 0
    
    A[0][1] = A[0][1] / A[0][0]
    B[0] = B[0] / A[0][0]
    for i in range (1, n - 1):
        A[i][i + 1] = A[i][i + 1] / (A[i][i] - A[i + 1][i]*A[i - 1][i])
        B[i] = (B[i] - B[i - 1] *  A[i + 1][i])/(A[i][i] -  A[i + 1][i]*A[i - 1][i])

    y_sec[n] = B[n - 1]
    for i in range (n - 1, 0, -1):
        y_sec[i] = B[i - 1] - A[i - 1][i]* y_sec[i + 1]
    return y_sec

def dichotomic_search (x, x_array, start, end):
    if(end - start <= 1):
        return start
    index = (start+end)//2
    if x < x_array[index] :
        return dichotomic_search(x, x_array, start, index)
    elif x > x_array[index]: 
        return dichotomic_search(x, x_array, index, end)
    else:
        return index

def zeta(x, x_array, i):
    return (x - x_array[i])/(x_array[i + 1] - x_array[i])

def interpolation_result(x, x_array, y_array, n, y_sec, i):
    if (i == n - 1):
        return y_array[n - 1]
    k = x_array[i + 1] - x_array[i]
    z = zeta(x, x_array, i)
    
    inter_x = (1 - z) * y_array[i]
    inter_x = inter_x + z * y_array[i + 1]
    inter_x = inter_x + ((k**2)/6) * ((1 - z)**3 - (1 - z)) * y_sec[i]
    inter_x = inter_x + ((k**2)/6) * (z**3 - z) * y_sec[i + 1]
    return inter_x

def interpolation(x_array, y_array):
    n = x_array.size
    A = create_matrix_A(x_array, n)
    B = create_matrix_B(x_array, y_array, n)
    y_sec = solve_linear_equations(A, B, n - 2)
    f = (lambda x: interpolation_result(x, x_array, y_array, n, y_sec, dichotomic_search(x, x_array, 0, n)))
    return f

    
    

    
