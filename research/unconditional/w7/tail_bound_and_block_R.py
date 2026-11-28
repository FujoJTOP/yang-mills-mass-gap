import numpy as np
from scipy.special import iv

beta = 2.5

# (1) Explicit tail bound for M_f:
#     (M_f)[a,b] = I_{|a-b|}(lam) - I_{a+b+2}(lam),  |(M_f)[a,b]| <= I_{|a-b|}(lam)
#     sum_{d>=D} I_d(lam) <= e^{lam^2/4} * (lam/2)^D / D! / (1 - lam/(2D))   for D > lam/2
def tail_bound(lam, D):
    if D <= lam / 2:
        return np.inf
    return np.exp(lam ** 2 / 4) * (lam / 2) ** D / np.math.factorial(D) / (1 - lam / (2 * D))

import math
def tail_bound2(lam, D):
    if D <= lam / 2:
        return np.inf
    return np.exp(lam ** 2 / 4) * (lam / 2) ** D / math.factorial(D) / (1 - lam / (2 * D))

print("=== explicit tail bound vs actual sum_{d>=D} I_d(lam) ===")
for lam in (1.0, 1.5, 3.0):
    for D in (6, 8, 10, 12):
        actual = sum(iv(d, lam) for d in range(D, 400))
        bnd = tail_bound2(lam, D)
        print("  lam=%.1f D=%2d : actual=%.3e  bound=%.3e  ratio=%.2f"
              % (lam, D, actual, bnd, bnd / actual if actual > 0 else np.inf))
    print()

print("=== relative error budget: 2*tail_bound / t_0(beta) with t_0 = 2 I_1(beta)/beta ===")
for beta in (1.0, 2.5, 4.0):
    t0 = 2 * iv(1, beta) / beta
    print("  beta=%.1f : t_0=%.6f" % (beta, t0))
    for lam in (1.0, 2.0, 3.0):
        for D in (8, 12, 16):
            b = 2 * tail_bound2(lam, D)
            print("     lam=%.1f D=%2d : 2*tail/t_0 = %.3e" % (lam, D, b / t0))
    print()

# (2) high-precision check of (R) on blocks
print("=== (R) on D x D blocks: r(T_D) vs I2/I1 ===")
def Tfull(beta, lam, J2):
    t = np.array([2 * iv(k + 1, beta) / beta for k in range(J2 + 1)])
    M = np.array([[iv(abs(a - b), lam) - iv(a + b + 2, lam) for b in range(J2 + 1)]
                  for a in range(J2 + 1)])
    return M @ np.diag(t) @ M

for beta in (0.5, 1.0, 2.5, 4.0, 6.0):
    rK = float(iv(2, beta) / iv(1, beta))
    print("  beta=%.1f : I2/I1 = %.10f" % (beta, rK))
    for lam in (0.25, 0.5, 1.0, 2.0, 3.0):
        for D in (6, 10, 14):
            T = Tfull(beta, lam, 40)
            ev = np.sort(np.linalg.eigvalsh((T[:D + 1, :D + 1] + T[:D + 1, :D + 1].T) / 2))[::-1]
            r = ev[1] / ev[0]
            print("     lam=%.2f D=%2d : r=%.10f  %s" % (lam, D, r, "OK" if r <= rK + 1e-12 else "FAIL"))
    print()
