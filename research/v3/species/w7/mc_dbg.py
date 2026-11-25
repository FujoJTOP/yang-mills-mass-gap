"""mc_dbg.py -- settle the insertion convention for one plaquette, numerically."""
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
    return Ts


LAM = [1j * T for T in su_generators(2)]


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


# word: A B C D  with D = (link_3)^{-1} etc.  Take all four links positive first.
A, B, C, D = haar(), haar(), haar(), haar()
M = [A, B, C, D]
h = 1e-3
for k in range(4):
    for a in range(3):
        W = list(M); W[k] = M[k] @ expL(h, a)
        fd = (np.trace(prod(W)).real - np.trace(prod(M)).real) / h
        rot = M[k:] + M[:k]
        head = np.trace(LAM[a] @ prod(rot)).real
        tail = np.trace(LAM[a] @ prod(rot[1:] + [rot[0]])).real
        print(f"k={k} a={a}  fd={fd:+.8f}  head={head:+.8f}  tail={tail:+.8f}")
    break
print()
for k in range(4):
    for a in range(3):
        W = list(M); W[k] = M[k] @ expL(h, a)
        fd = (np.trace(prod(W)).real - np.trace(prod(M)).real) / h
        rot = M[k:] + M[:k]
        head = np.trace(LAM[a] @ prod(rot)).real
        tail = np.trace(LAM[a] @ prod(rot[1:] + [rot[0]])).real
        print(f"k={k} a={a}  fd={fd:+.8f}  head={head:+.8f}  tail={tail:+.8f}")
