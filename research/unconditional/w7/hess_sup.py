import numpy as np

def make_QP(q, a, beta=1.0):
    a = np.asarray(a, float)
    def QP(Y):
        Y = np.asarray(Y, float).reshape(4, 3)
        s = Y.sum(axis=0)
        tot = beta * q * (s @ s)
        cr = np.zeros(3)
        for i in range(4):
            for j in range(i + 1, 4):
                cr += np.cross(Y[i], Y[j])
        return tot - 2.0 * beta * (cr @ a)
    return QP

def hess(QP, n=12):
    H = np.zeros((n, n))
    for k in range(n):
        e = np.zeros(n); e[k] = 1.0
        H[k, k] = QP(e)
    for k in range(n):
        for l in range(k + 1, n):
            e = np.zeros(n); e[k] = 1.0; e[l] = 1.0
            ek = np.zeros(n); ek[k] = 1.0
            el = np.zeros(n); el[l] = 1.0
            H[k, l] = H[l, k] = (QP(e) - QP(ek) - QP(el)) / 2.0
    return H

print("=== sup over (q,a), q^2+a^2=1, of |lambda_min| / beta ===")
best = 0.0; bestq = None
for th in np.linspace(0, np.pi / 2, 2001):
    q, a = np.cos(th), np.sin(th)
    H = hess(make_QP(q, [0, 0, a]))
    ev = np.linalg.eigvalsh(H)
    r = -ev.min()
    if r > best:
        best = r; bestq = (q, a, ev.min(), ev.max())
print("  sup = %.6f at q=%.4f, a=%.4f  (lambda_min=%.4f, lambda_max=%.4f)"
      % (best, bestq[0], bestq[1], bestq[2], bestq[3]))
print("  1+sqrt(2) = %.6f" % (1 + np.sqrt(2)))
print()
print("=== per-edge Hessian bound (6 plaquettes per edge) and the resulting beta range ===")
C = 6 * best
print("  ||Hess S|| <= 6 * %.6f * beta = %.6f * beta" % (best, C))
print("  H >= (2 - %.4f beta) I" % C)
print("  H >= I  <=>  beta <= 1/(%.4f - ?) ... solve 2 - %.4f*beta >= 1  =>  beta <= %.6f"
      % (C, C, 1.0 / C))
print("  (old bound 30 beta  =>  beta <= 1/30 = %.6f)" % (1 / 30))
print("  improvement factor = %.3f x" % ((1 / C) / (1 / 30)))
print()
print("=== check: gauge directions are exactly zero for all (q,a) ===")
worst = 0.0
for th in np.linspace(0, np.pi / 2, 51):
    q, a = np.cos(th), np.sin(th)
    QP = make_QP(q, [0, 0, a])
    for d in range(3):
        Yg = np.zeros(12)
        for i, s in enumerate([1, -1, 1, -1]):
            Yg[3 * i + d] = s
        worst = max(worst, abs(QP(Yg)))
print("  max |Q_P(gauge dir)| over all (q,a) and d = %.3e" % worst)
