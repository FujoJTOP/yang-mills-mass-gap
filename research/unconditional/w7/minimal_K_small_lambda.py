import numpy as np
from scipy.special import iv

# Round 40: close the small-lambda gap of the single-link (R) block route.
# For each lambda find the minimal K such that
#       r_K + eps_K / T00  <=  I2/I1
# If K(lambda) stays bounded as lambda -> 0, the route closes.

beta = 2.5
rK0 = float(iv(2, beta) / iv(1, beta))

def build(lam, J2=80):
    M = np.array([[iv(abs(a - b), lam) - iv(a + b + 2, lam) for b in range(J2 + 1)]
                  for a in range(J2 + 1)])
    t = np.array([2 * iv(k + 1, beta) / beta for k in range(J2 + 1)])
    T = M @ np.diag(t) @ M
    return (T + T.T) / 2

def try_K(T, K, T00):
    b = T[:K, :K]
    ev = np.sort(np.linalg.eigvalsh((b + b.T) / 2))[::-1]
    r_k = ev[1] / ev[0]
    P = np.zeros_like(T); P[:K, :K] = b
    eps = np.linalg.norm(T - P, 2)
    return r_k + eps / T00, r_k, eps

print("=== minimal K such that  r_K + eps_K/T00 <= I2/I1 ===")
print("  beta=%.1f  I2/I1=%.12f" % (beta, rK0))
print("%10s %10s %8s %16s %16s %14s" %
      ("lam", "I2/I1-r_K", "K_min", "r_K", "eps_K/T00", "lhs"))
for lam in (0.0, 1e-6, 1e-4, 1e-3, 1e-2, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0):
    T = build(lam)
    T00 = T[0, 0]
    found = None
    for K in range(4, 41):
        lhs, r_k, eps = try_K(T, K, T00)
        if lhs <= rK0:
            found = (K, r_k, eps, lhs); break
    if found is None:
        # report the best K
        best = min((try_K(T, K, T00)[0], K) for K in range(4, 41))
        print("%10.0e %10s %8s %16s %16s %14.10f  NO K<=40" %
              (lam, "-", "-", "-", "-", best[0]))
    else:
        K, r_k, eps, lhs = found
        print("%10.0e %10.6f %8d %16.10f %16.3e %14.10f" %
              (lam, rK0 - r_k, K, r_k, eps / T00, lhs))
