import numpy as np
from scipy.special import iv

# CANDIDATE EXACT IDENTITY discovered numerically:
#       M_f D M_f = D          equivalently   M_f^{-1} = D M_f D
#   with M_f[a,b] = I_{|a-b|}(lam) - I_{a+b+2}(lam),   D = diag((-1)^a).
# Verify and probe.

def Mf(J2, lam):
    return np.array([[iv(abs(a - b), lam) - iv(a + b + 2, lam) for b in range(J2 + 1)]
                     for a in range(J2 + 1)])

for lam in (0.3, 0.8, 1.5, 2.5, 4.0):
    for J2 in (6, 10, 16):
        M = Mf(J2, lam)
        D = np.diag([(-1.0) ** a for a in range(J2 + 1)])
        L = M @ D @ M
        res = np.abs(L - D).max()
        print("lam=%.1f J2=%2d : max |M_f D M_f - D| = %.3e" % (lam, J2, res))

print()
print("=== if the identity holds, then (M_f D)^2 = I, so spec(M_f D) = +-1 ===")
lam = 1.5; J2 = 12
M = Mf(J2, lam); D = np.diag([(-1.0) ** a for a in range(J2 + 1)])
w = np.linalg.eigvalsh(M @ D)
print("  eigenvalues of M_f D :", np.array2string(np.sort(w), precision=8))
print("  away from +-1 by at most %.3e" % np.min(np.abs(np.abs(w) - 1)))
print()
print("=== consequence: T^{-1} = M_f^{-1} K^{-1} M_f^{-1} = D M_f K^{-1} M_f D ===")
beta = 2.5
t = np.array([2 * iv(m + 1, beta) / beta for m in range(J2 + 1)])
K = np.diag(t)
T = M @ K @ M
Tinv = np.linalg.inv(T)
rhs = D @ M @ np.linalg.inv(K) @ M @ D
print("  max |T^{-1} - D M_f K^{-1} M_f D| = %.3e" % np.abs(Tinv - rhs).max())
