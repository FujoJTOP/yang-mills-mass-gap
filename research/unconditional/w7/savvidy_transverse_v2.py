import numpy as np
from scipy.linalg import expm

# (a) step 1, corrected:
#   S_P = beta (1 - (1/2) Re tr( U1 U2 U3^-1 U4^-1 ))
#   Savvidy point: U1 = i sigma3, U2 = U3 = U4 = I   (=> q_P = 0, |a_P| = 1)
#   12 real variables eps = (eps_1..eps_4) in R^3 each,  U_i = exp(i eps_i.sigma/2) U_i^0
#   Gauge orbit: U_e -> g_a U_e g_b^-1 with (g_0,g_1,g_2,g_3) = 12 params.
#   Project the 12x12 Hessian onto the orthogonal complement of the 12 gauge tangents,
#   and report the physical (transverse) spectrum.

s1 = np.array([[0,1],[1,0]], dtype=complex)
s2 = np.array([[0,-1j],[1j,0]], dtype=complex)
s3 = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2, dtype=complex)
sig = [s1, s2, s3]

def su2(v):
    return expm(1j * (v[0]*s1 + v[1]*s2 + v[2]*s3) / 2.0)

U0 = [1j*s3, I2.copy(), I2.copy(), I2.copy()]
beta = 1.0

def Sp_from_U(Us, beta=beta):
    P = Us[0] @ Us[1] @ np.linalg.inv(Us[2]) @ np.linalg.inv(Us[3])
    return beta * (1.0 - 0.5*np.real(np.trace(P)))

def Sp(eps, beta=beta):
    return Sp_from_U([su2(eps[3*i:3*i+3]) @ U0[i] for i in range(4)], beta)

def hessian(f, n=12, h=2e-4):
    H = np.zeros((n, n))
    f0 = f(np.zeros(n))
    for a in range(n):
        ea = np.zeros(n); ea[a] = h
        for b in range(a, n):
            eb = np.zeros(n); eb[b] = h
            eab = ea + eb
            val = (f(eab) - f(ea) - f(eb) + f0) / (h*h)
            H[a, b] = H[b, a] = val
    return H

H = hessian(Sp)
print("=== 12x12 Hessian at Savvidy (beta=1), h=2e-4 ===")
ev = np.sort(np.linalg.eigvalsh(H))
print("  ", np.array2string(ev, precision=6, max_line_width=220))

# gauge tangents: 4 g's acting as U_e -> g_a U_e g_b^-1
# link i connects corner a_i (start) and b_i (end): U1: 0->1, U2: 1->2, U3: 3->2, U4: 0->3
links = [(0,1), (1,2), (3,2), (0,3)]
hg = 1e-6
G = np.zeros((12,12))
col = 0
for corner in range(4):
    for k in range(3):
        phi = np.zeros(3); phi[k] = hg
        g = su2(phi); gi = np.linalg.inv(g)
        gs = [I2.copy() for _ in range(4)]
        gs[corner] = g
        gis = [I2.copy() for _ in range(4)]
        gis[corner] = gi
        for i, (a_i, b_i) in enumerate(links):
            Unew = gs[a_i] @ U0[i] @ gis[b_i]
            dU = Unew - U0[i]
            M = dU @ np.linalg.inv(U0[i])
            v = np.array([-1j*np.real(np.trace(M @ sig[a])) for a in range(3)])
            G[3*i:3*i+3, col] = v / hg
        col += 1

print()
print("=== gauge tangent matrix: 12 vectors; rank and orthonormality ===")
U_sv, S_sv, Vt = np.linalg.svd(G, full_matrices=False)
print("  singular values:", np.array2string(S_sv, precision=6))
rank = int((S_sv > 1e-8).sum())
print("  numerical rank = %d  (12 params, degenerate at this point => lowering of rank)" % rank)

Q = U_sv[:, :rank]
P = np.eye(12) - Q @ Q.T
Hp = P @ H @ P
evp = np.sort(np.linalg.eigvalsh(Hp))
print()
print("=== projected Hessian (gauge removed) : 12 eigenvalues ===")
print("  ", np.array2string(evp, precision=8, max_line_width=220))
dim_phys = 12 - rank
print()
print("=== physical (transverse) dimension = %d ; its spectrum ===" % dim_phys)
phys = evp[:dim_phys] if dim_phys > 0 else np.array([])
print("  ", np.array2string(phys, precision=8, max_line_width=220))
if dim_phys > 0:
    print("  min = %.10f" % phys.min())
