import numpy as np
from scipy.special import iv

# Round 46: which term carries the ordering  N_c(y2) <= N_c(y1) ?
#   N_c = (1/4)[ 2 + (W+ + W- - y0^2 t0 - E+ - E-)/D ]
#   normalise: wT = (W+ + W-)/D ,  yT = y0^2 t0 / D ,  eT = (E+ + E-)/D
#   => N_c = (1/4)[ 2 + wT - yT - eT ]

def build(beta, lam, J=140):
    a = np.arange(J + 1)
    M = iv(np.abs(a[:, None] - a[None, :]), lam) - iv(a[:, None] + a[None, :] + 2, lam)
    t = 2 * iv(a + 1, beta) / beta
    T = M @ np.diag(t) @ M
    return (T + T.T) / 2, M, t

def parts(y, t):
    Wp = np.sum(y[:-1] ** 2 * t[1:])
    Wm = np.sum(y[1:] ** 2 * t[:-1])
    Ep = np.sum((y[1:] - y[:-1]) ** 2 * t[1:])
    Em = np.sum((y[1:] - y[:-1]) ** 2 * t[:-1])
    D = np.sum(y ** 2 * t)
    yT = y[0] ** 2 * t[0] / D
    wT = (Wp + Wm) / D
    eT = (Ep + Em) / D
    nc = 0.25 * (2 + wT - yT - eT)
    return nc, wT, yT, eT

beta = 2.5
print("beta = %.1f" % beta)
print()
print("N_c = (1/4)[2 + wT - yT - eT] ;  ordering needs (wT-yT-eT) smaller for y2")
print()
print("%6s | %10s %10s %10s %10s | %10s %10s %10s %10s" %
      ("lam", "Nc(y1)", "wT1", "yT1", "eT1", "Nc(y2)", "wT2", "yT2", "eT2"))
for lam in (0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
    T, M, t = build(beta, lam)
    ev, V = np.linalg.eigh(T)
    y1, y2 = M @ V[:, -1], M @ V[:, -2]
    n1, w1, yo1, e1 = parts(y1, t)
    n2, w2, yo2, e2 = parts(y2, t)
    print("%6.2f | %10.6f %10.6f %10.6f %10.6f | %10.6f %10.6f %10.6f %10.6f" %
          (lam, n1, w1, yo1, e1, n2, w2, yo2, e2))
print()
print("=== which pieces separate?  (need wT2-wT1 <= (yT2-yT1)+(eT2-eT1) ... i.e. eT does the work) ===")
print("%6s %14s %14s %14s %14s" % ("lam", "wT2-wT1", "yT2-yT1", "eT2-eT1", "LHS-RHS"))
for lam in (0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
    T, M, t = build(beta, lam)
    ev, V = np.linalg.eigh(T)
    y1, y2 = M @ V[:, -1], M @ V[:, -2]
    _, w1, yo1, e1 = parts(y1, t)
    _, w2, yo2, e2 = parts(y2, t)
    lhs = (w2 - w1) - (yo2 - yo1)
    rhs = e2 - e1
    print("%6.2f %14.6f %14.6f %14.6f %14.6f %s" % (lam, w2 - w1, yo2 - yo1, e2 - e1, lhs - rhs,
                                                    "OK" if lhs - rhs <= 0 else "FAIL"))
