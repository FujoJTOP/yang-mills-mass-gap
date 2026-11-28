import numpy as np
from scipy.special import iv
from scipy.linalg import expm

# EXACT SU(2) single-link model:
#   character basis chi_j,  j = 0, 1/2, 1, ..., jmax   (dim = 2*jmax+1)
#   K is diagonal (convolution):  K_{jj} = t_j = 2*I_{2j+1}(beta)/beta
#   q = (1/2) tr U = (1/2) chi_{1/2};  q * chi_j = (1/2)(chi_{j-1/2} + chi_{j+1/2})
#       => M_q is tridiagonal with 1/2 off-diagonals
#   f = exp(-beta'*(1-q))  <= 1   =>  M_f = exp(-beta'*(I - M_q))
#   T = M_f K M_f
# Check claim (c):  <phi_1, f^2 phi_1> >= <phi_2, f^2 phi_2>

def model(beta, betap, jmax=8):
    js = np.arange(0, jmax + 1) / 2.0          # j = 0, .5, 1, ...
    d = len(js)
    t = np.array([2.0 * iv(2 * j + 1, beta) / beta for j in js])   # K eigenvalues
    K = np.diag(t)
    # M_q tridiagonal
    Mq = np.zeros((d, d))
    for i in range(d):
        if i - 1 >= 0: Mq[i, i - 1] = 0.5
        if i + 1 < d:  Mq[i, i + 1] = 0.5
    # M_f = exp(-betap*(I - Mq))
    Mf = expm(-betap * (np.eye(d) - Mq))
    T = Mf @ K @ Mf
    w, V = np.linalg.eigh(T)
    idx = np.argsort(-w); w = w[idx]; V = V[:, idx]
    p1, p2 = V[:, 0], V[:, 1]
    def w2(p):
        fp = Mf @ p
        return fp @ fp
    return w2(p1), w2(p2), w[0], w[1], t[1] / t[0]

print("=== EXACT SU(2) single-link model: check (c)  ||f phi_2|| <= ||f phi_1|| ===")
print("%7s %7s %7s %14s %14s %10s %10s %10s %6s"
      % ("beta", "beta'", "jmax", "||f p1||^2", "||f p2||^2", "n2/n1", "r(K)", "r(T)", "(c)"))
bad = 0; tot = 0
for jmax in (6, 10, 14):
    for beta in (0.5, 1.0, 2.5, 4.0):
        for betap in (0.5, 1.0, 2.5, 5.0):
            try:
                n1, n2, t1, t2, rK = model(beta, betap, jmax)
            except Exception as e:
                continue
            ratio = n2 / n1 if n1 > 0 else np.inf
            ok = ratio <= 1 + 1e-9
            tot += 1
            if not ok: bad += 1
            rT = t2 / t1 if t1 > 0 else np.nan
            print("%7.2f %7.2f %7d %14.8f %14.8f %10.6f %10.6f %10.6f %6s"
                  % (beta, betap, jmax, n1, n2, ratio, rK, rT, "OK" if ok else "FAIL"))
    print()
print("  (c)-failures = %d / %d" % (bad, tot))
