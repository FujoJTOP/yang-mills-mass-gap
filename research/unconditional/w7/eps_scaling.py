import numpy as np
from scipy.special import iv

# Round 43: identify the true controlling quantity of eps_K.
# Suspects:  t_K/t_0  and  I_K(lambda).

beta = 40.0
lam = 1.0
J2 = 200

def build(beta, lam, J2):
    M = np.array([[iv(abs(a - b), lam) - iv(a + b + 2, lam) for b in range(J2 + 1)]
                  for a in range(J2 + 1)])
    t = np.array([2 * iv(k + 1, beta) / beta for k in range(J2 + 1)])
    T = M @ np.diag(t) @ M
    return (T + T.T) / 2, M, t

T, M, t = build(beta, lam, J2)
t0 = t[0]
T00 = T[0, 0]
print("beta=%.1f lam=%.1f  t_0=%.6e  T00=%.6e" % (beta, lam, t0, T00))
print()
print("%5s %14s %14s %14s %14s %14s" %
      ("K", "eps_K/T00", "t_K/t_0", "I_K(lam)", "I_K^2", "eps/t_K"))
for K in (8, 12, 16, 24, 32, 48, 64):
    P = np.zeros_like(T); P[:K, :K] = T[:K, :K]
    eps = np.linalg.norm(T - P, 2)
    print("%5d %14.6e %14.6e %14.6e %14.6e %14.6e"
          % (K, eps / T00, t[K] / t0, iv(K, lam), iv(K, lam) ** 2, eps / t[K]))

print()
print("=== where does the corner mass actually sit? (T entries by row index) ===")
print("%5s %14s %14s %14s" % ("a", "T[a,a]/T00", "max_b |T[a,b]|/T00", "sum_b T[a,b]^2/T00^2"))
for a in (0, 4, 8, 12, 16, 24, 32, 48):
    if a < T.shape[0]:
        print("%5d %14.6e %14.6e %14.6e"
              % (a, T[a, a] / T00, np.abs(T[a, :]).max() / T00, np.sum(T[a, :] ** 2) / T00 ** 2))
