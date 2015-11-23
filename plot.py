from matplotlib import pyplot as plt
import numpy as np
import iterative
import pascal
import power
plt.style.use('ggplot')
qr = []
lu = []

for i in range(2, 13):
    q = pascal.solve_qr_b(pascal.pascal_matrix(i), pascal.harmonic_vector(i))
    l = pascal.solve_lu_b(pascal.pascal_matrix(i), pascal.harmonic_vector(i))
    qr.append(q)
    lu.append(l)

plt.subplot(1, 1, 1)
x = range(2, 13)
y = [i[1] for i in qr]
z = [i[2] for i in qr]

plt.plot(x, y, color='blue')  # error from householder
plt.plot(x, z, color='green')  # solution error of qr
plt.yscale('log')
plt.savefig('./qr_err.png')

y = [i[1] for i in lu]
z = [i[2] for i in lu]
plt.clf()
plt.plot(x, y, color='blue')
plt.plot(x, z, color='green')
plt.yscale('log')
plt.savefig('./lu_err.png')
plt.clf()


jacobi, gs = iterative.generate_data()
j_vals = [i[1] for i in jacobi]
g_vals = [i[1] for i in gs]

jacobi_approx = sum(j_vals) / len(j_vals)  # 2c
gs_approx = sum(g_vals) / len(g_vals)

print("Averages, jacobi then gauss-seidel, then iterations")
print(jacobi_approx)
print(gs_approx)
print(float(sum(j[2] for j in jacobi))/sum(g[2] for g in gs))

exact = np.array([9.0/190, 28.0/475, 33.0/475]).reshape(3,1)
errs_jacobi = [pascal.norm_inf(j-exact) for j in j_vals]
errs_gs = [pascal.norm_inf(g-exact) for g in g_vals]

plt.plot([j[2] for j in jacobi], errs_jacobi, 'ko', [g[2] for g in gs], errs_jacobi, 'bo')
plt.savefig('./iterative_err.png')
plt.clf()


powers = power.generate_data()
ds = [p[0] for p in powers if p[0] is not None]
ts = [p[1] for p in powers if p[1] is not None]
tis = [p[2] for p in powers if p[2] is not None]
maxs = [p[3] for p in powers if p[3] is not None]
mins = [p[4] for p in powers if p[4] is not None]
big = max(maxs)
small = max(mins)
maxs = [float(m)/big for m in maxs]
mins = [float(m)/small for m in mins]


plt.scatter(ds, ts, c=maxs)
plt.savefig('./power_mat.png')
plt.clf()
plt.scatter([1.0/d for d in ds], tis, c=mins)
plt.savefig('./power_inv.png')
plt.clf()