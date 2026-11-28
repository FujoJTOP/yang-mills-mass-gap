import numpy as np
from scipy.special import iv

# Reduction idea: T = M_f K M_f.  How well is T approximated by its top-left KxK block?
# If ||T - P_K T P_K|| is small compared with t_1, then (R) reduces to a FINITE
# computation (KxK block with explicit Bessel entries).

beta = 2.5
J2 = 60
t = np.array([2 * iv(k + 1, beta) / beta for k in range(J2 + 1)])

def build(lam, J2=J2):
    M = np.array([[iv(abs(a - b), lam) - iv(a + b + 2, lam) for b in range(J2 + 1)]
                  for a in range(J2 + 1)])
    return M @ np.diag(t) @ M

print("=== ||T - P_K T P_K|| / t_1  and the ratio error ===")
print("%6s %10s %14s %14s %12s" % ("K", "lam", "||T-PTP||/t1", "r(true)", "r(K-block)"))
for lam in (0.5, 1.5, 3.0, 6.0):
    T = build(lam)
    ev = np.sort(np.linalg.eigvalsh((T + T.T) / 2))[::-1]
    rtrue = ev[1] / ev[0]
    for K in (2, 4, 6, 8, 12, 16):
        P = np.zeros_like(T); P[:K, :K] = T[:K, :K]
        err = np.linalg.norm(T - P, 2) / ev[0]
        evb = np.sort(np.linalg.eigvalsh((P[:K, :K] + P[:K, :K].T) / 2))[::-1]
        print("%6d %10.1f %14.3e %14.8f %12.8f   err_r=%+.2e"
              % (K, lam, err, rtrue, evb[1] / evb[0], evb[1] / evb[0] - rtrue))
    print()

print("=== how fast do the entries of M_f decay away from the diagonal (bandwidth)? ===")
lam = 1.5
M = np.array([[iv(abs(a - b), lambda_ := lam) - iv(a + b + 2, lam) for b in range(41)] for a in range(41)])
print("  row 20, |M_f| by |a-b|:")
for d in (0, 1, 2, 4, 8, 12, 16, 20):
    print("    |a-b|=%2d : %.6e  (relative %.2e)" % (d, abs(M[20, 20 + d]) if 20 + d < 41 else 0,
                                                     abs(M[20, 20 + d]) / M[20, 20] if 20 + d < 41 else 0))
print()
print("=== same for T ===")
T = build(lam, 20)
print("  row 10, |T| by |a-b| (relative to T[10,10]=%.6e):" % T[10, 10])
for d in (0, 1, 2, 3, 5, 8, 10):
    v = abs(T[10, 10 + d]) if 10 + d < 21 else 0
    print("    |a-b|=%2d : %.6e  (relative %.3e)" % (d, v, v / T[10, 10]))
