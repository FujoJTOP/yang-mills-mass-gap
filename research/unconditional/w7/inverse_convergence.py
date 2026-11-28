import numpy as np
from scipy.special import iv

# Does M_f^{-1} converge to the pure Toeplitz( I_{|a-b|} ) as J2 -> infinity?

def Mf(J2, lam):
    M = np.zeros((J2 + 1, J2 + 1))
    for a in range(J2 + 1):
        for b in range(J2 + 1):
            M[a, b] = iv(abs(a - b), lam) - iv(a + b + 2, lam)
    return M

lam = 1.5
print("=== convergence of (M_f^{-1})_{0,0} and (M_f^{-1})_{0,1} ===")
print("   I_0(1.5) = %.12f ;  I_1(1.5) = %.12f" % (iv(0, lam), iv(1, lam)))
for J2 in (6, 10, 16, 24, 34, 46):
    M = Mf(J2, lam)
    Mi = np.linalg.inv(M)
    print("  J2=%2d : (M^{-1})_{00} = %.12f   (M^{-1})_{01} = %.12f   (M^{-1})_{0,6} = %.6e"
          % (J2, Mi[0, 0], Mi[0, 1], Mi[0, min(6, J2)]))

print()
print("=== alternatively: compute M_f^{-1} for the BILATERAL Toeplitz part ===")
print("   check the Neumann-series identity  sum_{m in Z} I_m(x) I_{m+k}(x) = I_k(2x) ?")
x = 1.5
for k in (0, 1, 2):
    s = sum(iv(m, x) * iv(m + k, x) for m in range(-300, 301))
    print("   k=%d : sum = %.12f   I_k(2x) = %.12f" % (k, s, iv(k, 2 * x)))

print()
print("=== and: is M_f^{1/2} helpful? compute M_f^{-1/2} bands ===")
M = Mf(16, lam)
w, V = np.linalg.eigh(M)
Mih = V @ np.diag(1.0 / np.sqrt(np.maximum(w, 1e-300))) @ V.T
for k in range(0, 6):
    band = np.array([Mih[i, i + k] for i in range(17 - k)])
    print("  |i-j|=%d : rel %.6f (approx I_%d(1.5)/I_0 = %.6f)"
          % (k, band.max() / Mih[0, 0], k, iv(k, lam) / iv(0, lam)))
