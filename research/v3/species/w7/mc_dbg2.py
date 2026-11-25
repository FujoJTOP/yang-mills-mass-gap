"""mc_dbg2.py -- diagonal block of D2 for a word with inverted factors: which insertion order?"""
import numpy as np
rng = np.random.default_rng(7)
N = 2


def su_generators(N):
    Ts = []
    for j in range(N):
        for k in range(j + 1, N):
            T = np.zeros((N, N), dtype=complex); T[j, k] = T[k, j] = 0.5; Ts.append(T)
    for j in range(N):
        for k in range(j + 1, N):
            T = np.zeros((N, N), dtype=complex); T[j, k] = -0.5j; T[k, j] = 0.5j; Ts.append(T)
    for l in range(1, N):
        T = np.zeros((N, N), dtype=complex)
        d = 1.0 / np.sqrt(2.0 * l * (l + 1))
        for j in range(l):
            T[j, j] = d
        T[l, l] = -l * d
        Ts.append(T)
    return Ts


LAM = [1j * T for T in su_generators(2)]
NG = len(LAM)
print("NG =", NG)


def haar():
    Z = (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) / np.sqrt(2.0)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R)
    Q = Q * (d / np.abs(d)).conj()[None, :]
    return Q / np.linalg.det(Q) ** (1.0 / N)


def expL(t, a):
    return np.cos(t / 2) * np.eye(2) + 2.0 * np.sin(t / 2) * LAM[a]


def prod(l):
    P = np.eye(2, dtype=complex)
    for m in l:
        P = P @ m
    return P


U = [haar() for _ in range(4)]
signs = [+1, -1, +1, -1]
h = 1e-3


def word(Us, sg):
    return [Us[i] if sg[i] > 0 else Us[i].conj().T for i in range(4)]


for k in range(4):
    rot = word(U, signs)[k:] + word(U, signs)[:k]
    head = prod(rot)
    tail = prod(rot[1:] + [rot[0]])
    for a in range(NG):
        for b in range(NG):
            def val(sa, sb):
                V = list(U); V[k] = U[k] @ expL(sa, a) @ expL(sb, b)
                return np.trace(prod(word(V, signs))).real
            fd = (val(h, h) - val(h, 0) - val(0, h) + val(0, 0)) / h ** 2
            c_ab_tail = np.trace(LAM[a] @ LAM[b] @ tail).real
            c_ba_tail = np.trace(LAM[b] @ LAM[a] @ tail).real
            c_ab_head = np.trace(LAM[a] @ LAM[b] @ head).real
            c_ba_head = np.trace(LAM[b] @ LAM[a] @ head).real
            best = min([('ab_tail', c_ab_tail), ('ba_tail', c_ba_tail),
                        ('ab_head', c_ab_head), ('ba_head', c_ba_head)],
                       key=lambda t: abs(t[1] - fd))
            print(f"  sign={signs[k]:+d} k={k} a={a} b={b} fd={fd:+9.6f} "
                  f"ab_tail={c_ab_tail:+9.6f} ba_tail={c_ba_tail:+9.6f} "
                  f"ab_head={c_ab_head:+9.6f} ba_head={c_ba_head:+9.6f}   best={best[0]}")
        if k > 0:
            break
    if k > 1:
        break
