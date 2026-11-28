import numpy as np
from scipy.special import iv

# (1) How much of the top two eigenvectors of T sits in the span of the first few characters?
# (2) Is r(T) well approximated by the top-left 2x2 block ratio?

beta = 2.5
J2 = 20

def Mf(J2, lam):
    return np.array([[iv(abs(a - b), lam) - iv(a + b + 2, lam) for b in range(J2 + 1)]
                     for a in range(J2 + 1)])

def build(beta, lam, J2=J2):
    t = np.array([2 * iv(m + 1, beta) / beta for m in range(J2 + 1)])   # t_{j=m/2} = 2 I_{m+1}/beta
    M = Mf(J2, lam)
    return M @ np.diag(t) @ M, t

print("=== eigenvector concentration on the first few characters ===")
for lam in (0.0, 0.5, 1.5, 3.0):
    T, t = build(beta, lam)
    ev, V = np.linalg.eigh((T + T.T) / 2)
    idx = np.argsort(-ev); ev = ev[idx]; V = V[:, idx]
    r = ev[1] / ev[0]
    print("  lam=%.1f : r=%.8f" % (lam, r))
    for k in (1, 2, 4, 8):
        m1 = np.sum(V[:k, 0] ** 2); m2 = np.sum(V[:k, 1] ** 2)
        print("     first %2d comps: mass(x1)=%.6f  mass(x2)=%.6f" % (k, m1, m2))
    # 2x2 block prediction
    Tb = T[:2, :2]
    evb = np.sort(np.linalg.eigvalsh((Tb + Tb.T) / 2))[::-1]
    print("     2x2 block ratio = %.8f   (vs true %.8f)" % (evb[1] / evb[0], r))
    # 3x3 and 4x4
    for k in (3, 4, 6):
        Tb = T[:k, :k]
        evb = np.sort(np.linalg.eigvalsh((Tb + Tb.T) / 2))[::-1]
        print("     %dx%d block ratio = %.8f" % (k, k, evb[1] / evb[0]))
    print()

print("=== the 2x2 block explicitly ===")
for lam in (0.0, 0.5, 1.5, 3.0):
    T, t = build(beta, lam)
    print("  lam=%.1f  T[:2,:2] = [[%.8f, %.8f],[%.8f, %.8f]]  t[:2]=[%.8f,%.8f]"
          % (lam, T[0, 0], T[0, 1], T[1, 0], T[1, 1], t[0], t[1]))
