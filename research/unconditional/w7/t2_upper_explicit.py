import numpy as np
from scipy.special import iv

# Single-link model:  T = M_f K M_f  in the character basis (j = k/2).
#   (M_f)[a,b] = I_{|a-b|}(lam) - I_{a+b+2}(lam)      (a,b = doubled indices)
#   t_j = 2 I_{j2+1}(beta)/beta ,  K = diag
#
# Route for (R) = [t2/t1 <= I2/I1]:
#   (i)  t1 >= T00  where T00 = <u_0, T u_0> = T[0,0]      -- PROVED (round 360)
#   (ii) need an UPPER bound on t2 of the form  t2 <= (I2/I1) * T00.
#   Explicit test object: take g = the top eigenvector of the leading KxK block of T
#   (embedded), then  t2 <= max_{x perp g} R(x) = ||P_g T P_g||_op .
# Test whether (ii) holds with such g.

def build(beta, lam, J2=40):
    M = np.array([[iv(abs(a - b), lam) - iv(a + b + 2, lam) for b in range(J2 + 1)]
                  for a in range(J2 + 1)])
    t = np.array([2 * iv(k + 1, beta) / beta for k in range(J2 + 1)])
    T = M @ np.diag(t) @ M
    return (T + T.T) / 2, t, M

def top2(T):
    ev = np.sort(np.linalg.eigvalsh(T))[::-1]
    return ev[0], ev[1]

print("=== single-link model: can an explicit g give t2 <= (I2/I1)*T00 ? ===")
print("%6s %6s %14s %14s %14s %14s %10s" %
      ("beta", "lam", "T00", "t1", "t2", "t2_ub(g)", "(I2/I1)T00"))
for beta in (2.5, 4.0):
    rK = float(iv(2, beta) / iv(1, beta))
    for lam in (0.3, 0.6, 1.0, 1.5, 2.5):
        T, t, M = build(beta, lam, J2=40)
        t1, t2 = top2(T)
        T00 = T[0, 0]
        best = None
        for K in (2, 3, 4, 6, 8):
            ev, V = np.linalg.eigh(T[:K, :K])
            g = np.zeros(T.shape[0]); g[:K] = V[:, -1]
            g = g / np.linalg.norm(g)
            P = np.eye(T.shape[0]) - np.outer(g, g)
            Tp = P @ T @ P
            ub = np.sort(np.linalg.eigvalsh((Tp + Tp.T) / 2))[-1]
            if best is None or ub < best[1]:
                best = (K, ub)
        print("%6.1f %6.1f %14.8f %14.8f %14.8f %14.8f %10.8f   K=%d  ub<=rK*T00? %s"
              % (beta, lam, T00, t1, t2, best[1], rK * T00, best[0],
                 "YES" if best[1] <= rK * T00 + 1e-12 else "NO"))
