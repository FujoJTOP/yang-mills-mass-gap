import numpy as np
from scipy.special import iv
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
import functools, time

# FIX (round 58): in the earlier two-link runs the CHARACTER SUM was truncated to J <= jmax2
# at the same time as the BASIS.  These are different truncations:
#   - the basis restriction j <= Jout  is the (legitimate) compression PMP,
#   - the character sum must be carried to convergence (J <= Jchar ~ 6 at betap=2).
# This script separates the two.

@functools.lru_cache(maxsize=None)
def d3(j1, m1, j2, m2, j3, m3):
    return float(sN(wigner_3j(Rational(j1, 2), Rational(j2, 2), Rational(j3, 2),
                             Rational(m1, 2), Rational(m2, 2), Rational(m3, 2)), 25))

def make_basis(jmax2):
    return [(j2, m2, n2) for j2 in range(0, jmax2 + 1)
            for m2 in range(-j2, j2 + 1, 2) for n2 in range(-j2, j2 + 1, 2)]

def mult_op(B, index, J2, M2, N2, conj=False):
    d = len(B); A = np.zeros((d, d))
    mm, nn = (-M2, -N2) if conj else (M2, N2)
    ph = ((-1.0) ** ((M2 - N2) // 2)) if conj else 1.0
    for i2, (j2, m2, n2) in enumerate(B):
        m3, n3 = m2 - mm, n2 - nn
        for j3 in range(abs(J2 - j2), J2 + j2 + 1, 2):
            key = (j3, m3, n3)
            if key not in index: continue
            val = ((-1.0) ** ((m2 - n2) // 2)) * d3(j3, m3, J2, mm, j2, -m2) * d3(j3, n3, J2, nn, j2, -n2)
            if val != 0.0:
                A[index[key], i2] = np.sqrt((j3 + 1.0) * (j2 + 1.0)) * val
    return ph * A

def run(beta, Jout, Jchar):
    B = make_basis(Jout); d = len(B); index = {b: i for i, b in enumerate(B)}
    n = d * d
    bp = beta
    Mf = np.zeros((n, n))
    t0 = time.time()
    for J2 in range(0, Jchar + 1):
        cJ = (1.0 if J2 == 0 else 2 * (J2 + 1) * iv(J2 + 1, bp) / bp)
        for M2 in range(-J2, J2 + 1, 2):
            for N2 in range(-J2, J2 + 1, 2):
                Mf += cJ * np.kron(mult_op(B, index, J2, M2, N2, False),
                                   mult_op(B, index, J2, M2, N2, True))
    Mf *= np.exp(-bp)
    tj = np.array([2 * iv(j2 + 1, beta) / beta for j2 in range(Jout + 1)])
    Kd = np.array([tj[b1[0]] * tj[b2[0]] for b1 in B for b2 in B])
    T = (Mf * Kd[None, :]) @ Mf
    T = (T + T.T) / 2
    ev = np.sort(np.linalg.eigvalsh(T))[::-1]
    return ev[1] / ev[0], ev[0], time.time() - t0

rK = float(iv(2, 2.5) / iv(1, 2.5))
print("=== two-link, beta=2.5 : separate the character truncation from the basis ===")
print("%6s %8s %16s %16s %10s" % ("Jout", "Jchar", "r(T)", "lambda1", "sec"))
for Jout in (2, 3):
    for Jchar in (Jout, 4, 6):
        r, l1, dt = run(2.5, Jout, Jchar)
        print("%6d %8d %16.10f %16.6e %10.1f" % (Jout, Jchar, r, l1, dt))
    print()
print("  I2/I1 = %.10f" % rK)
