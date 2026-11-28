import numpy as np
from scipy.special import iv
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
import functools, sys

# Optimized exact two-link SU(2), j <= 3/2 (dim 30 per link, 900 total).
# Goal: verify whether r(T) really turns over at large beta (round 363 needed j<=1 only).

@functools.lru_cache(maxsize=None)
def d3(j1, m1, j2, m2, j3, m3):
    return float(sN(wigner_3j(Rational(j1, 2), Rational(j2, 2), Rational(j3, 2),
                             Rational(m1, 2), Rational(m2, 2), Rational(m3, 2)), 25))

def make_basis(jmax2):
    B = []
    for j2 in range(0, jmax2 + 1):
        for m2 in range(-j2, j2 + 1, 2):
            for n2 in range(-j2, j2 + 1, 2):
                B.append((j2, m2, n2))
    return B

def mult_op_fast(B, index, J2, M2, N2, conj=False):
    d = len(B)
    A = np.zeros((d, d))
    mm, nn = (-M2, -N2) if conj else (M2, N2)
    ph = ((-1.0) ** ((M2 - N2) // 2)) if conj else 1.0
    for i2, (j2, m2, n2) in enumerate(B):
        m3, n3 = m2 - mm, n2 - nn
        if m3 < 0 or n3 < 0:
            continue
        lo = abs(J2 - j2)
        for j3 in range(lo, J2 + j2 + 1, 2):
            key = (j3, m3, n3)
            if key not in index:
                continue
            i3 = index[key]
            sgn = (-1.0) ** ((m2 - n2) // 2)
            val = sgn * d3(j3, m3, J2, mm, j2, -m2) * d3(j3, n3, J2, nn, j2, -n2)
            if val == 0.0:
                continue
            A[i3, i2] = np.sqrt((j3 + 1.0) * (j2 + 1.0)) * val
    return ph * A

def main(jmax2, betas, betap_fixed=None):
    B = make_basis(jmax2)
    d = len(B)
    index = {b: i for i, b in enumerate(B)}
    print("jmax2=%d : dim per link = %d ; two links = %d" % (jmax2, d, d * d), flush=True)
    for beta in betas:
        bp = beta if betap_fixed is None else betap_fixed
        Mf = np.zeros((d * d, d * d))
        scale = np.exp(-bp)
        for J2 in range(0, jmax2 + 1):
            cJ = (1.0 if J2 == 0 else (2 * (J2 + 1) * iv(J2 + 1, bp) / bp if bp > 0 else 0.0))
            for M2 in range(-J2, J2 + 1, 2):
                for N2 in range(-J2, J2 + 1, 2):
                    A1 = mult_op_fast(B, index, J2, M2, N2, False)
                    A2 = mult_op_fast(B, index, J2, M2, N2, True)
                    Mf += cJ * np.kron(A1, A2)
        Mf *= scale
        tj = {j2: 2 * iv(j2 + 1, beta) / beta for j2 in range(jmax2 + 1)}
        K = np.diag([tj[b1[0]] * tj[b2[0]] for b1 in B for b2 in B])
        T = Mf @ K @ Mf
        T = (T + T.T) / 2
        ev = np.sort(np.linalg.eigvalsh(T))[::-1]
        r = ev[1] / ev[0]
        print("  beta=%6.1f betap=%6.1f : r=%.10f  lam1=%.6e  1-r=%.6f  beta(1-r)=%.4f"
              % (beta, bp, r, ev[0], 1 - r, beta * (1 - r)), flush=True)
        del Mf, K, T

print("=== sanity first: betap=0 must give M_f = I and r = I2/I1 ===")
main(3, [2.5], betap_fixed=0.0)
print()
print("=== physical line betap = beta, j <= 3/2 ===")
main(3, [0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0])

