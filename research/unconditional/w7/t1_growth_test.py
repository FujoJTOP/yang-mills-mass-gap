import numpy as np
from scipy.special import iv

# REDUCTION found in round 20:
#   T(lambda) = M_w K_op M_w  <=  K_op  (w <= 1 pointwise), hence
#       t_2(T) <= t_2(K_op)  automatically.
#   Therefore  (R)  follows from the SINGLE inequality
#       t_1(T(lambda)) >= t_1(K_op) = t_0 = 2 I_1(beta)/beta .
# Verify numerically over a range of beta and lambda.

beta = 2.5
N = 3000
th = (np.arange(N) + 0.5) * np.pi / N
dth = np.pi / N
mu = (2 / np.pi) * dth

def K_closed(a, b, beta):
    return (np.exp(beta * np.cos(a - b)) - np.exp(beta * np.cos(a + b))) / (2 * beta)

TH, PH = np.meshgrid(th, th, indexing='ij')

def top_eigs(beta, lam, k=3):
    Kg = K_closed(TH, PH, beta)
    w = np.exp(lam * np.cos(th))
    A = Kg * mu
    T = np.diag(w) @ A @ np.diag(w)
    wv = np.sort(np.linalg.eigvalsh((T + T.T) / 2))[::-1]
    return wv[:k], np.sort(np.linalg.eigvalsh((A + A.T) / 2))[::-1][:k]

print("=== t_1(T(lambda)) vs t_1(K_op) = 2 I_1(beta)/beta ===")
for beta in (1.0, 2.5, 4.0):
    t0 = 2 * iv(1, beta) / beta
    print("  beta=%.1f : t_1(K_op) = %.10f" % (beta, t0))
    for lam in (0.0, 0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 4.0):
        tT, tK = top_eigs(beta, lam)
        ok = tT[0] >= tK[0] - 1e-9
        print("     lam=%.2f : t_1(T)=%.8f  t_1(K)=%.8f  diff=%+.3e  %s"
              % (lam, tT[0], tK[0], tT[0] - tK[0], "OK" if ok else "FAIL"))
    print()

print("=== and the resulting ratio bound ===")
for beta in (1.0, 2.5, 4.0):
    t0 = 2 * iv(1, beta) / beta
    t1K = 2 * iv(2, beta) / beta
    print("  beta=%.1f : I2/I1 = %.10f" % (beta, t1K / t0))
    for lam in (0.5, 1.5, 3.0):
        tT, tK = top_eigs(beta, lam)
        print("     lam=%.1f : t_2(T)/t_1(T) = %.8f  vs I2/I1 = %.8f  %s"
              % (lam, tT[1] / tT[0], t1K / t0, "OK" if tT[1] / tT[0] <= t1K / t0 + 1e-9 else "FAIL"))
