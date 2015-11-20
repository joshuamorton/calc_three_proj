import pascal
import numpy as np


def power_method(A, v, epsilon, max_iter):
    """
    return lambda, e or (None, None)
    """
    u_0 = v
    # u_0 = u_0 / pascal.norm(u_0)
    w = np.array([1, 1]).reshape(2, 1)

    iteration = 0
    u_n = u_0
    while iteration < max_iter:
        u_1 = pascal.mult(A, u_n)
        # u_1 = u_1 / pascal.norm(u_1)  # normalize
        error = pascal.norm_inf(u_1 / pascal.norm(u_1) - u_n / pascal.norm(u_n))
        if epsilon > error:
            l = pascal.mult(w.T, u_1) / pascal.mult(w.T, u_n) 
            u_1 = u_1 / pascal.norm(u_1)
            u_n = u_n / pascal.norm(u_n)
            return u_1, l, iteration
        u_n = u_1
        iteration += 1
    else:
        return None, None, None


def rand_matrix():
    return np.random.rand(2, 2) * 4 - 2


def inv(A):
    a, b, c, d = A[0, 0], A[0, 1], A[1, 0], A[1, 1]
    if (a*d - b * c) == 0:
        return None
    det = 1/(a*d - b*c)

    return np.matrix([[d, -b], [-c, a]]) * det, det


def generate_data():
    i = 0
    data_array = []
    while i < 1000:
        mat = rand_matrix()
        mat_inv, determinant = inv(mat)
        if mat_inv is None:
            continue
        u_max, l_max, iterations_max = power_method(mat,
                                        np.array([1, 1]).reshape(2, 1), .00005,
                                        100)

        u_min, l_min, iterations_min = power_method(mat,
                                        np.array([1, 1]).reshape(2, 1), .00005,
                                        100)
        trace = np.trace(mat)
        if all([determinant, trace, iterations_max, iterations_min]):
            data_array.append((determinant, trace, iterations_max, iterations_min))
            i+= 1
    return data_array

