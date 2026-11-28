import numpy as np
from scipy.special import iv
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN

# EXACT two-link SU(2) Peter-Weyl with the ORTHONORMAL basis e_{jmn} = sqrt(2j+1) D^j_{mn}.
#   (M_f)_{i3,i2} = sqrt((2j3+1)(2j2+1)) * int D_{i3} f conj(D_{i2}) dU
#   K = K1 (x) K1 diagonal with t_{j1} t_{j2},  t_j = 2 I_{2j+1}(beta)/beta
#   T = M_f K M_f ;  r(T) = t2/t1
# Sanity: betap = 0 => f = 1 => M_f = I.

def d3(j1, m1, j2, m2, j3, m3):
    return float(sN(wigner_3j(Rational(j1,2), Rational(j2,2), Rational(j3,2),
                             Rational(m1,2), Rational(m2,2), Rational(m3,2)), 30))

def make_basis(jmax2):
    B = []
    for j2 in range(0, jmax2 + 1):
        for m2 in range(-j2, j2 + 1, 2):
            for n2 in range(-j2, j2 + 1, 2):
                B.append((j2, m2, n2))
    return B

def mult_op_orth(B, J2, M2, N2, conj=False):
    """orthonormal-basis matrix of multiplication by D^J_{MN} (or conj), restricted to B."""
    d = len(B)
    A = np.zeros((d, d))
    mm, nn = (-M2, -N2) if conj else (M2, N2)
    ph = ((-1.0) ** ((M2 - N2) // 2)) if conj else 1.0
    for i3, (j3, m3, n3) in enumerate(B):
        for i2, (j2, m2, n2) in enumerate(B):
            if m3 != m2 - mm or n3 != n2 - nn:
                continue
            sgn = (-1.0) ** ((m2 - n2) // 2)
            val = sgn * d3(j3, m3, J2, mm, j2, -m2) * d3(j3, n3, J2, nn, j2, -n2)
            if val == 0.0:
                continue
            A[i3, i2] = np.sqrt((j3 + 1.0) * (j2 + 1.0)) * val   # sqrt((2j+1)(2j'+1))
    return ph * A

def build(beta, betap, jmax2=2):
    B = make_basis(jmax2)
    d = len(B)
    tj = {j2: 2 * iv(j2 + 1, beta) / beta for j2 in range(jmax2 + 1)}
    Mf = np.zeros((d * d, d * d))
    scale = np.exp(-betap)
    for J2 in range(0, jmax2 + 1):
        cJ = (1.0 if J2 == 0 else 2 * (J2 + 1) * iv(J2 + 1, betap) / betap) if betap > 0 else (1.0 if J2 == 0 else 0.0)
        for M2 in range(-J2, J2 + 1, 2):
            for N2 in range(-J2, J2 + 1, 2):
                A1 = mult_op_orth(B, J2, M2, N2, conj=False)
                A2 = mult_op_orth(B, J2, M2, N2, conj=True)
                Mf += cJ * np.kron(A1, A2)
    Mf *= scale
    K = np.diag([tj[b1[0]] * tj[b2[0]] for b1 in B for b2 in B])
    T = Mf @ K @ Mf
    return T, Mf, B, tj

print("=== sanity: betap=0 => M_f = I ? ===")
T, Mf, B, tj = build(2.5, 0.0, jmax2=2)
print("  dim=%d ; max |M_f - I| = %.3e" % (Mf.shape[0], np.abs(Mf - np.eye(Mf.shape[0])).max()))

print()
print("=== sanity: M_f symmetric ?  and T at betap=0 should be K ===")
T0, Mf0, B, tj = build(2.5, 0.0, jmax2=2)
K0 = np.diag([tj[b1[0]] * tj[b2[0]] for b1 in B for b2 in B])
print("  max |T - K| at betap=0 = %.3e" % np.abs(T0 - K0).max())

print()
print("=== two-link ratio ===")
for betap in (0.0, 0.3, 0.8, 1.5, 3.0):
    T, Mf, B, tj = build(2.5, betap, jmax2=2)
    T = (T + T.T) / 2
    ev = np.sort(np.linalg.eigvalsh(T))[::-1]
    rK = float(iv(2, 2.5) / iv(1, 2.5))
    print("  betap=%.1f : top4 = %s" % (betap, np.array2string(ev[:4], precision=8)))
    print("      r(T)=%.10f  I2/I1=%.10f  %s"
          % (ev[1] / ev[0], rK, "OK" if ev[1] / ev[0] <= rK + 1e-9 else "FAIL"))



