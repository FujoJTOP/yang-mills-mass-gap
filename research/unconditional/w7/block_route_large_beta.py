import numpy as np
from scipy.special import iv

# Round 41: does the block route survive large beta?
# eps_K decays in K at a rate controlled by t_K = 2 I_{K+1}(beta)/beta, and for n << beta
# Bessel functions are NOT small. So K may need to grow with beta.

def build(beta, lam, J2=100):
    M = np.array([[iv(abs(a - b), lam) - iv(a + b + 2, lam) for b in range(J2 + 1)]
                  for a in range(J2 + 1)])
    t = np.array([2 * iv(k + 1, beta) / beta for k in range(J2 + 1)])
    T = M @ np.diag(t) @ M
    return (T + T.T) / 2

def lhs_for(T, K, T00):
    b = T[:K, :K]
    ev = np.sort(np.linalg.eigvalsh((b + b.T) / 2))[::-1]
    r_k = ev[1] / ev[0]
    P = np.zeros_like(T); P[:K, :K] = b
    eps = np.linalg.norm(T - P, 2)
    return r_k + eps / T00, r_k, eps

print("=== does r_K + eps_K/T00 <= I2/I1 survive large beta ? ===")
for beta in (2.5, 5.0, 10.0, 20.0, 40.0):
    rK0 = float(iv(2, beta) / iv(1, beta))
    lam = 1.0
    T = build(beta, lam, J2=120)
    T00 = T[0, 0]
    print("  beta=%5.1f lam=%.1f  I2/I1=%.10f  T00=%.4e" % (beta, lam, rK0, T00))
    for K in (8, 12, 16, 24, 32, 48, 64):
        lhs, r_k, eps = lhs_for(T, K, T00)
        print("     K=%3d : r_K=%.10f  eps/T00=%.3e  lhs=%.10f  %s"
              % (K, r_k, eps / T00, lhs, "OK" if lhs <= rK0 + 1e-12 else "FAIL"))
    print()
