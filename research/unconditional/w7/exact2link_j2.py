import numpy as np
from scipy.special import iv
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
import functools, gc, time

# (C) CONVERGENCE TEST: exact two-link SU(2), j <= 2 (dim 55 per link, 3025 total).
# Question left open by rounds 362/364: is the truncation converged? j<=1 gave 0.569 at beta=2,
# j<=3/2 gave 0.880. Now j<=2.

@functools.lru_cache(maxsize=None)
def d3(j1, m1, j2, m2, j3, m3):
    return float(sN(wigner_3j(Rational(j1, 2), Rational(j2, 2), Rational(j3, 2),
                             Rational(m1, 2), Rational(m2, 2), Rational(m3, 2)), 25))

def make_basis(jmax2):
    return [(j2, m2, n2)
            for j2 in range(0, jmax2 + 1)
            for m2 in range(-j2, j2 + 1, 2)
            for n2 in range(-j2, j2 + 1, 2)]

def mult_op_fast(B, index, J2, M2, N2, conj=False):
    d = len(B); A = np.zeros((d, d))
    mm, nn = (-M2, -N2) if conj else (M2, N2)
    ph = ((-1.0) ** ((M2 - N2) // 2)) if conj else 1.0
    for i2, (j2, m2, n2) in enumerate(B):
        m3, n3 = m2 - mm, n2 - nn
        lo = abs(J2 - j2)
        for j3 in range(lo, J2 + j2 + 1, 2):
            key = (j3, m3, n3)
            if key not in index: continue
            i3 = index[key]
            val = ((-1.0) ** ((m2 - n2) // 2)) * d3(j3, m3, J2, mm, j2, -m2) * d3(j3, n3, J2, nn, j2, -n2)
            if val != 0.0:
                A[i3, i2] = np.sqrt((j3 + 1.0) * (j2 + 1.0)) * val
    return ph * A

def run(jmax2, betas):
    B = make_basis(jmax2); d = len(B); n = d * d
    index = {b: i for i, b in enumerate(B)}
    print("jmax2=%d : d=%d, two-link dim=%d, M_f is %.1f MB" % (jmax2, d, n, n * n * 8 / 1e6), flush=True)
    for beta in betas:
        t0 = time.time()
        bp = beta
        Mf = np.zeros((n, n))
        for J2 in range(0, jmax2 + 1):
            cJ = (1.0 if J2 == 0 else 2 * (J2 + 1) * iv(J2 + 1, bp) / bp)
            for M2 in range(-J2, J2 + 1, 2):
                for N2 in range(-J2, J2 + 1, 2):
                    A1 = mult_op_fast(B, index, J2, M2, N2, False)
                    A2 = mult_op_fast(B, index, J2, M2, N2, True)
                    Mf += cJ * np.kron(A1, A2)
        Mf *= np.exp(-bp)
        tj = np.array([2 * iv(j2 + 1, beta) / beta for j2 in range(jmax2 + 1)])
        Kd = np.array([tj[b1[0]] * tj[b2[0]] for b1 in B for b2 in B])
        T = Mf * Kd[None, :]           # = Mf @ diag(Kd)
        T = T @ Mf
        T = (T + T.T) / 2
        ev = np.sort(np.linalg.eigvalsh(T))[::-1]
        r = ev[1] / ev[0]
        print("  beta=%5.1f betap=%5.1f : r=%.10f  I2/I1=%.10f  1-r=%.6f  [%.0fs]"
              % (beta, bp, r, float(iv(2, beta) / iv(1, beta)), 1 - r, time.time() - t0), flush=True)
        del Mf, T; gc.collect()

print("=== j<=2 means jmax2=4 : d=55, two-link dim=3025 ===")
run(4, [2.0])

