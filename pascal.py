import numpy as np
import math


def combination(n, r):
    return math.factorial(n) / math.factorial(r) / math.factorial(n - r)


def pascal_matrix(n):
    """
    computes a pascal matrix such as
    [[1,1,1,1],
     [1,2,3,4],
     [1,3,6,10],
     [1,4,10,20]]
    """
    return np.matrix([[combination(x+y-2, y-1) for y in range(1,n+1)]
                      for x in range(1,n+1)], dtype='double')

def harmonic_vector(n):
    """
    create a vector in the form [1,1/2,1/3,...1/n]
    """
    return np.array([[1.0 / i] for i in range(1, n + 1)], dtype='double')


def mult(a, b):
    """
    assumes the given matrices are the correct shape
    """
    result = np.zeros((a.shape[0], b.shape[1]))
    for i in range(result.shape[0]):
        for j in range(result.shape[1]):
            # result[i,j] = a col i * b row j
            result[i,j] = np.dot(a[i], b[:,j])
    return result


def norm_inf(matrix):
    return max(abs(sum(np.nditer(row))) for row in matrix)

def norm(vector):
    return math.sqrt(sum(float(x)**2 for x in np.nditer(vector)))


def lu_fact(matrix):  # 1a
    """
    computes the LU factorization of a matrix
    """
    size = len(matrix)
    l_mat = np.zeros((size,size))
    u_mat = np.zeros((size,size)) + matrix
    for row in range(size):
        for col in range(size):
            if row == col:
                l_mat[row,col] = 1 / u_mat[row,col]
                u_mat[row] /= u_mat[row,col]
            elif row > col:
                l_mat[row, col] += u_mat[row,col] / u_mat[col,col]
                u_mat[row] -= (u_mat[row,col] / u_mat[col,col]) * u_mat[col]
    
    err = norm_inf(mult(l_mat, u_mat) - matrix)

    return l_mat, u_mat, err


def qr_fact_househ(matrix):  # 1b
    size = matrix.shape[1]
    big = max(matrix.shape)
    working = np.zeros(matrix.shape) + matrix  # copy
    hs = []
    for col in range(size-1):
        submatrix = working[col:,col:]  # work only on the submatrix
        v = submatrix[:,0]  # gets a column
        e = np.zeros(submatrix.shape[1])  # calculates e for the submatrix
        e[0] = 1
        u = v + norm(v) * e
        u /= norm(u)  # calculate u = (v + ||v|| * e) / ||(v + ||v|| * e)||
        H_part = np.eye(submatrix.shape[0]) - 2 * mult(np.matrix(u).T, np.matrix(u))
        H = np.zeros((big,big))
        for i in range(col):
            H[i,i] = 1
        H[col:,col:] = H_part
        working = mult(H, working)
        hs.append(H)

    q_mat = reduce(mult, hs, np.eye(size))
    r_mat = mult(reduce(mult, hs[::-1], np.eye(size)), matrix)
    err = norm_inf(mult(q_mat, r_mat) - matrix)

    return q_mat, r_mat, err


def qr_fact_givens(matrix):  # 1b
    size = matrix.shape[1]
    working = np.zeros(matrix.shape) + matrix
    gs = []
    for col in range(size):
        for row in range(col+1, size):
            x, y = working[col, col], working[row, col]
            sqrtx2y2 = math.sqrt(x**2 + y**2)
            c = x / sqrtx2y2
            s = -y / sqrtx2y2
            givensmat = np.eye(size)
            givensmat[col, col] = c
            givensmat[row, row] = c
            givensmat[col, row] = -s
            givensmat[row, col] = s
            working = mult(givensmat, working)
            gs.append(givensmat)

    q = reduce(mult, (g.T for g in gs), np.eye(size))
    r = mult(reduce(mult, gs[::-1], np.eye(size)), matrix)
    err = norm_inf(mult(q, r) - matrix)
    return q, r, err


def solve_lu_b(A, b):  # 1c
    l, u, err = lu_fact(A)
    b_copy = b.copy()
    c_vec = [0 for n in b]
    for row in range(len(b)):
        curr = b_copy[row]
        for col in range(row):
            curr -= l[row, col] * c_vec[col]
        c_vec[row] = curr / l[row, row]

    x_vec = [0 for n in b]
    for row in range(len(b)-1, -1, -1):
        curr = c_vec[row]
        col = len(b) - 1
        while col > row:
            curr -= u[row, col] * x_vec[col]
            col -= 1
        x_vec[row] = curr / u[row, row]

    x_vec = np.array(x_vec)
    sol_err = norm_inf(mult(A, x_vec) - b)
    return x_vec, err, sol_err


def solve_qr_b(A, b):  # 1c
    q, r, err = qr_fact_househ(A)
    c_vec = mult(q.T, b)

    x_vec = [0 for n in b]
    for row in range(len(b)-1, -1, -1):
        curr = c_vec[row]
        col = len(b) - 1
        while col > row:
            curr -= r[row, col] * x_vec[col]
            col -= 1
        x_vec[row] = curr / r[row, row]

    x_vec = np.array(x_vec)
    sol_err = norm_inf(mult(A, x_vec) - b)
    return x_vec, err, sol_err


def load_matrix(file):
    with open(file) as f:
        return np.matrix([line.split(",") for line in f.readlines()], dtype='double')


def convert_matrix(m):
    A = m[:, :-1]
    b = m[:, -1]
    return A, b


def main():  # 1d
    print "solving with lu"
    for i in range(2, 13):
        x, err, sol_err = solve_lu_b(pascal_matrix(i), harmonic_vector(i))
        print i
        print x
        print "norm(LU-P), norm(Px-b)"
        print err, sol_err
        print ""


    print "solving with qr"
    for i in range(2, 13):
        x, err, sol_err = solve_qr_b(pascal_matrix(i), harmonic_vector(i))
        print i
        print x
        print "norm(QR-P), norm(Px-b)"
        print err, sol_err
        print ""

