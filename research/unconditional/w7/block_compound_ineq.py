import numpy as np
from scipy.special import iv

# Single-link (R) via block truncation:
#   t1(T) >= T00                        [T], round 360
#   r(T) <= r_K + eps_K / t1(T)         r_K = t2/t1 of the KxK leading block,
#                                        eps_K = ||T - P_K T P_K||_op
#   t1(T) >= T00  =>  eps_K/t1(T) <= eps_K/T00
#   => (R) follows from   r_K + eps_K/T00 <= I2/I1
# Test this compound inequality over lambda.

def build(beta, lam, J2=60):
    M = np.array([[iv(abs(a - b), lam) - iv(a + b + 2, lam) for b in range(J2 + 1)]
                  for a in range(J2 + 1)])
    t = np.array([2 * iv(k + 1, beta) / beta for k in range(J2 + 1)])
    T = M @ np.diag(t) @ M
    return (T + T.T) / 2

def block_stats(T, K):
    b = T[:K, :K]
    ev = np.sort(np.linalg.eigvalsh((b + b.T) / 2))[::-1]
    P = np.zeros_like(T); P[:K, :K] = b
    eps = np.linalg.norm(T - P, 2)
    return ev[0], ev[1], eps

print("=== compound inequality  r_K + eps_K/T00 <= I2/I1  (single-link model) ===")
beta = 2.5
rK = float(iv(2, beta) / iv(1, beta))
print("  beta=%.1f  I2/I1 = %.10f" % (beta, rK))
print("%8s %4s %14s %14s %14s %14s %10s" %
      ("lam", "K", "T00", "r_K", "eps_K/T00", "lhs", "verdict"))
for lam in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 1.5, 2.0, 2.5, 3.5, 5.0):
    T = build(beta, lam)
    T00 = T[0, 0]
    for K in (8, 12, 16):
        t1b, t2b, eps = block_stats(T, K)
        r_k = t2b / t1b
        lhs = r_k + eps / T00
        print("%8.2f %4d %14.8f %14.8f %14.3e %14.8f %10s"
              % (lam, K, T00, r_k, eps / T00, lhs, "OK" if lhs <= rK + 1e-12 else "FAIL"))
    print()

print("=== control: at lam = 0 the inequality must be an EQUALITY ===")
T = build(beta, 0.0)
print("  T00=%.10f  r_K (K=8)=%.10f  eps=%.3e  I2/I1=%.10f"
      % (T[0, 0], block_stats(T, 8)[1] / block_stats(T, 8)[0], block_stats(T, 8)[2], rK))
