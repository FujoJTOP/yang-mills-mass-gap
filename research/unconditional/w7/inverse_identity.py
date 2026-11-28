import numpy as np
from scipy.special import iv

# Apparent exact identity from the numerical inverse:
#       M_f^{-1} = Toeplitz( I_{|a-b|}(lambda) )     ??
# i.e.   M_f . T0 = I   with  M_f[a,b] = I_{|a-b|} - I_{a+b+2},  T0[b,c] = I_{|b-c|}
# Verify directly.

def Mf(J2, lam):
    M = np.zeros((J2 + 1, J2 + 1))
    for a in range(J2 + 1):
        for b in range(J2 + 1):
            M[a, b] = iv(abs(a - b), lam) - iv(a + b + 2, lam)
    return M

def T0(J2, lam):
    T = np.zeros((J2 + 1, J2 + 1))
    for a in range(J2 + 1):
        for b in range(J2 + 1):
            T[a, b] = iv(abs(a - b), lam)
    return T

for J2 in (6, 8, 10):
    lam = 1.5
    M = Mf(J2, lam); T = T0(J2, lam)
    prod = M @ T
    off = prod - np.eye(J2 + 1)
    print("J2=%2d lam=%.1f : max |M_f T0 - I| = %.3e   (max |M_f T0| = %.3e)"
          % (J2, lam, np.abs(off).max(), np.abs(prod).max()))

print()
print("=== diagonal / off-diagonal of M T0 ===")
for J2 in (10, 14):
    lam = 1.5
    M = Mf(J2, lam); T = T0(J2, lam)
    prod = M @ T
    diag = np.diag(prod)
    print("  J2=%2d : diag min=%.10f max=%.10f" % (J2, diag.min(), diag.max()))
    print("        offdiag max = %.3e" % np.abs(prod - np.diag(diag)).max())

print()
print("=== selected entries (J2=10, lam=1.5) ===")
lam = 1.5; J2 = 10
M = Mf(J2, lam); T = T0(J2, lam); P = M @ T
for (i, j) in [(0, 0), (0, 1), (2, 3), (5, 5), (9, 9), (1, 8)]:
    print("  (M T0)[%d,%d] = %.10f" % (i, j, P[i, j]))
