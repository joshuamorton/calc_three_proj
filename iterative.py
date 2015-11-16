import numpy as np
import pascal


def jacobi_iter(x_0, epsilon, max_iter):
    A = np.matrix([[1.0, 1.0/2, 1.0/3], [1.0/2, 1.0, 1.0/4], [1.0/3, 1.0/4, 1.0]])
    b = np.array([.1, .1, .1]).reshape(3, 1)
    iteration = 0
    error = 1000 
    x_n = x_0
    while iteration < max_iter:
        d = np.diag(1 / np.diag(A))
        alu = A - np.diag(np.diag(A))
        x_1 = pascal.mult(d, b) + pascal.mult(pascal.mult(d, -alu), x_n)
        error = pascal.norm_inf(x_1 - x_n)
        x_n = x_1
        if epsilon > error:
            return x_n, iteration
        iteration += 1
    else:
        return None, None


def gs_iter(x_0, epsilon, max_iter):
    A = np.matrix([[1.0, 1.0/2, 1.0/3], [1.0/2, 1.0, 1.0/4], [1.0/3, 1.0/4, 1.0]])
    B = np.array([.1, .1, .1]).reshape(3, 1)
    S_inv = np.zeros((3, 3))
    a, b, c, d, e, f = A[0, 0], A[1, 0], A[1, 1], A[2, 0], A[2, 1], A[2, 2]
    S_inv[0, 0] = 1.0/a
    S_inv[1, 0] = -b/(a*c)
    S_inv[1, 1] = 1.0/c
    S_inv[2, 0] = (-c*d+b*e)/(a*c*f)
    S_inv[2, 1] = -e/(c*f)
    S_inv[2, 2] = 1.0/f

    S = np.tril(A)
    U = A - S

    iteration = 0
    error = 1000 
    x_n = x_0
    while iteration < max_iter:
        x_1 = pascal.mult(pascal.mult(S_inv, -U), x_n)
        x_1 += pascal.mult(S_inv, B)
        error = pascal.norm_inf(x_1 - x_n)
        x_n = x_1
        if epsilon > error:
            return x_n, iteration
        iteration += 1
    else:
        return None, None


def rand_vec():
    return np.random.rand(3, 1) * 2 - 1
