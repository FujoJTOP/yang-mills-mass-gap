import numpy as np
from scipy.special import iv

# Phase A1: verify the Abel decomposition of N_c and locate which term dominates.
#
#  T(lam) = M_w K M_w ,  K = diag(t_a),  t_a = 2 I_{a+1}(beta)/beta   (doubled index a=0,1,2,...)
#  M_c = (1/2) S_1  (tridiagonal, off-diagonal 1/2)
#  N_c(y) = <y, M_c K y> / <y, K y>
#         = [ (1/2) sum_a y_a ( t_{a+1} y_{a+1} + t_{a-1} y_{a-1} ) ] / sum_a y_a^2 t_a
#
#  Abel:  sum_a y_a y_{a+1} t_{a+1}
#       = (1/2) sum_a (y_a^2 + y_{a+1}^2) t_{a+1}  -  (1/2) sum_a (y_{a+1} - y_a)^2 t_{a+1}
#       =: W_plus - (1/2) E
#
#  so   2*<y,M_c K y> = (W_plus + W_minus) - (1/2)(E_plus + E_minus)
#  where W_plus = sum_a y_a^2 t_{a+1} ,  W_minus = sum_{a>=1} y_a^2 t_{a-1}
#        E_plus = sum_a (y_{a+1}-y_a)^2 t_{a+1} ,  E_minus = sum_{a>=1} (y_a - y_{a-1})^2 t_{a-1}

def build(beta, lam, J=140):
    a = np.arange(J + 1)
    M = iv(np.abs(a[:, None] - a[None, :]), lam) - iv(a[:, None] + a[None, :] + 2, lam)
    t = 2 * iv(a + 1, beta) / beta
    T = M @ np.diag(t) @ M
    return (T + T.T) / 2, M, t

def Nc_abel(y, t):
    """N_c via the Abel decomposition.

    <y, M_c K y> = (1/4)[ W+ + W- + 2D - y_0^2 t_0 - E+ - E- ]
    with D = sum y_a^2 t_a.  Verified on y = e_0, which must give 0.
    """
    W_plus = np.sum(y[:-1] ** 2 * t[1:])          # sum_{a=0}^{J-1} y_a^2 t_{a+1}
    W_minus = np.sum(y[1:] ** 2 * t[:-1])         # sum_{a=1}^{J} y_a^2 t_{a-1}
    E_plus = np.sum((y[1:] - y[:-1]) ** 2 * t[1:])
    E_minus = np.sum((y[1:] - y[:-1]) ** 2 * t[:-1])
    D = np.sum(y ** 2 * t)
    num = 0.25 * (W_plus + W_minus + 2 * D - y[0] ** 2 * t[0] - E_plus - E_minus)
    return num / D, (W_plus + W_minus) / (4 * D), (E_plus + E_minus) / (4 * D), D

def Nc_direct(y, Mc, t):
    Ky = t * y
    return (y @ (Mc @ Ky)) / (y @ Ky)

def Mc_mat(J):
    M = np.zeros((J + 1, J + 1))
    for a in range(J):
        M[a, a + 1] = 0.5
        M[a + 1, a] = 0.5
    return M

beta = 2.5
print("beta = %.1f" % beta)
print()
print("=== (1) Abel decomposition reproduces N_c ===")
print("%6s %18s %18s %18s %18s" % ("lam", "N_c direct", "N_c Abel", "W-term", "E-term"))
for lam in (0.05, 0.2, 0.5, 1.0, 2.0, 3.0):
    T, M, t = build(beta, lam)
    Mc = Mc_mat(T.shape[0] - 1)
    ev, V = np.linalg.eigh(T)
    x1, x2 = V[:, -1], V[:, -2]
    y1, y2 = M @ x1, M @ x2
    n1d = Nc_direct(y1, Mc, t)
    n1a, W1, E1, _ = Nc_abel(y1, t)
    print("%6.2f %18.10f %18.10f %18.10f %18.10f" % (lam, n1d, n1a, W1, E1))
print()
print("=== (2) the two terms for y1 and y2 (which dominates the ordering?) ===")
print("%6s %14s %14s %14s %14s %14s %14s" %
      ("lam", "N_c(y1)", "N_c(y2)", "W(y1)", "W(y2)", "E(y1)", "E(y2)"))
for lam in (0.05, 0.2, 0.5, 1.0, 2.0, 3.0):
    T, M, t = build(beta, lam)
    ev, V = np.linalg.eigh(T)
    y1, y2 = M @ V[:, -1], M @ V[:, -2]
    n1, W1, E1, d1 = Nc_abel(y1, t)
    n2, W2, E2, d2 = Nc_abel(y2, t)
    print("%6.2f %14.8f %14.8f %14.8f %14.8f %14.3e %14.3e" % (lam, n1, n2, W1, W2, E1, E2))
