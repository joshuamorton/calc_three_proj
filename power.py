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
