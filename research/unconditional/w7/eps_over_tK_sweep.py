import numpy as np
from scipy.special import iv

# Round 43 cont.: test the candidate explicit bound   eps_K <= 10 * t_K .
# If it holds, then (R) follows from   r_K + 10 t_K / t_0 <= I2/I1 ,
# which has huge margin everywhere.

def eps_over_tK(beta, lam, K, J2=200):
    M = np.array([[iv(abs(a - b), lam) - iv(a + b + 2, lam) for b in range(J2 + 1)]
                  for a in range(J2 + 1)])
    t = np.array([2 * iv(k + 1, beta) / beta for k in range(J2 + 1)])
    T = M @ np.diag(t) @ M
    T = (T + T.T) / 2
    P = np.zeros_like(T); P[:K, :K] = T[:K, :K]
    eps = np.linalg.norm(T - P, 2)
    return eps / t[K], T, t

print("=== eps_K / t_K  ===")
print("%6s %6s %6s %14s" % ("beta", "lam", "K", "eps_K/t_K"))
worst = 0.0
for beta in (1.0, 2.5, 5.0, 10.0, 20.0, 40.0):
    for lam in (0.3, 1.0, 2.0):
        for K in (8, 12, 16, 24, 32):
            v, T, t = eps_over_tK(beta, lam, K)
            worst = max(worst, v)
            flag = "" if v <= 10 else "   <-- exceeds 10"
            if K in (12, 24) and lam in (0.3, 1.0):
                print("%6.1f %6.1f %6d %14.6e%s" % (beta, lam, K, v, flag))
print()
print("  WORST eps_K/t_K over the whole sweep = %.6f  (need <= 10)" % worst)
