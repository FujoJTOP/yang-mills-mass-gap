import numpy as np

# exact Hessian of Q_P at a configuration, using Q_P(e_k), Q_P(e_k+e_l)
# Q_P(Y) = beta*q*|sum y_i|^2 - 2*beta*sum_{i<j}(y_i x y_j).a     (exact quadratic)

def make_QP(q, a, beta):
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
    base = np.zeros(n)
    for k in range(n):
        e = np.zeros(n); e[k] = 1.0
        H[k, k] = QP(e) - 2 * QP(base) + QP(base)   # = QP(e)
    for k in range(n):
        for l in range(k + 1, n):
            e = np.zeros(n); e[k] = 1.0; e[l] = 1.0
            ek = np.zeros(n); ek[k] = 1.0
            el = np.zeros(n); el[l] = 1.0
            # Y^T H Y for Y=e: H_kk; for Y=e_k+e_l: H_kk+H_ll+2H_kl
            v = (QP(e) - QP(ek) - QP(el)) / 2.0
            H[k, l] = H[l, k] = v
    return H

beta = 2.5
cases = [
    ("Savvidy  (q=0, |a|=1)",      0.0, [0, 0, 1.0]),
    ("worst   (|q|=4/5,|a|=3/5)",  0.8, [0, 0, 0.6]),
    ("ordered (q=1, a=0)       ",  1.0, [0, 0, 0.0]),
    ("random  (q=0.3,a=0.954)  ",  0.3, [0, 0, np.sqrt(1 - 0.09)]),
]
for name, q, a in cases:
    QP = make_QP(q, a, beta)
    H = hess(QP)
    ev = np.linalg.eigvalsh(H)
    print("=== %s ===" % name)
    print("   spectrum:", np.array2string(ev, precision=4, suppress_small=True))
    print("   min = %+9.6f   max = %+9.6f   (beta=%.2f)" % (ev.min(), ev.max(), beta))
    # gauge directions: (+v,-v,+v,-v) in link space
    print("   gauge modes Q_P(+v,-v,+v,-v) for v=e_d :", end=" ")
    for d in range(3):
        Yg = np.zeros(12)
        for i, s in enumerate([1, -1, 1, -1]):
            Yg[3 * i + d] = s
        print("%+.4f" % (QP(Yg) / (Yg ** 2).sum()), end="  ")
    print()
    print()
