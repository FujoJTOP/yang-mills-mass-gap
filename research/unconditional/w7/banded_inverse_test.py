import numpy as np
from scipy.special import iv

# Structural test: is M_f^{-1} (and hence T^{-1}) BANDED?
#   M[x,y] = sum_{t=|x-y|}^{x+y} c_t   (c_t = I_{2t} - I_{2t+2}, positive, decreasing)
# Such "Toeplitz-Hankel" matrices are sometimes inverses of banded matrices.
# If M^{-1} is tridiagonal then T^{-1} = M^{-1} K^{-1} M^{-1} is pentadiagonal,
# which would reduce the spectral problem enormously.

def Mmat(J2, lam):
    M = np.zeros((J2 + 1, J2 + 1))
    for ia in range(J2 + 1):
        for ib in range(J2 + 1):
            lo, hi = abs(ia - ib), ia + ib
            M[ia, ib] = sum(iv(m, lam) - iv(m + 2, lam) for m in range(lo, hi + 1, 2))
    return M

J2 = 8
lam = 1.5
M = Mmat(J2, lam)
Minv = np.linalg.inv(M)

print("=== M_f for lam=%.1f (J2=%d) ===" % (lam, J2))
print(np.array2string(M, precision=6, max_line_width=200))
print()
print("=== M_f^{-1} : look at the decay away from the diagonal ===")
D = np.abs(Minv)
scale = np.abs(np.diag(Minv)).max()
for k in range(0, J2 + 1):
    band = np.array([D[i, i + k] for i in range(J2 + 1 - k)])
    print("  |i-j|=%d : max |M^{-1}| = %.6e   (relative %.2e)" % (k, band.max(), band.max() / scale))

print()
print("=== is M^{-1} exactly tridiagonal up to truncation? (relative threshold sweep) ===")
for tol in (1e-2, 1e-4, 1e-6, 1e-8):
    off = 0
    for i in range(J2 + 1):
        for j in range(J2 + 1):
            if abs(i - j) > 1 and D[i, j] / scale > tol:
                off += 1
    print("  entries with |i-j|>1 and relative size > %.0e : %d" % (tol, off))

print()
print("=== T = M K M with K=diag(t_j); is T^{-1} banded? ===")
beta = 2.5
t = np.array([2 * iv(m + 1, beta) / beta for m in range(J2 + 1)])   # t_j with m=2j
T = M @ np.diag(t) @ M
Tinv = np.linalg.inv(T)
DT = np.abs(Tinv); sc = np.abs(np.diag(Tinv)).max()
for k in range(0, J2 + 1):
    band = np.array([DT[i, i + k] for i in range(J2 + 1 - k)])
    print("  |i-j|=%d : max = %.6e  (rel %.2e)" % (k, band.max(), band.max() / sc))
