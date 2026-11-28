import numpy as np
from scipy.special import iv
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN

# FIXED exact two-link SU(2) Peter-Weyl computation.
#   M_f = exp(-betap) sum_J c_J sum_{M,N} M1[D^J_{MN}] (x) M2[conj D^J_{MN}]
#   conj(D^J_{MN}) = (-1)^{M-N} D^J_{-M,-N}
#   => M2[conj D^J_{MN}] = (-1)^{M-N} M[D^J_{-M,-N}]      <-- the FIX
# Sanity: at betap = 0, f == 1 so M_f == I.

def d3(j1, m1, j2, m2, j3, m3):
    return float(sN(wigner_3j(Rational(j1), Rational(m1), Rational(j2), Rational(m2),
                             Rational(j3), Rational(m3)), 30))

def make_basis(jmax2):
    B = []
    for j2 in range(0, jmax2 + 1):
        for m2 in range(-j2, j2 + 1, 2):
            for n2 in range(-j2, j2 + 1, 2):
                B.append((j2, m2, n2))
    return B

def mult_op(B, J2, M2, N2, conj=False):
    """matrix of multiplication by D^J_{M N} (or its conjugate) on L^2(SU(2))."""
    d = len(B)
    A = np.zeros((d, d))
    mm, nn = (-M2, -N2) if conj else (M2, N2)
    ph = ((-1) ** ((M2 - N2) // 2)) if conj else 1.0
    for i3, (j3, m3, n3) in enumerate(B):
        for i2, (j2, m2, n2) in enumerate(B):
            if m3 != mm + (-m2) or n3 != nn + (-n2):
                continue
            sgn = (-1) ** ((m2 - n2) // 2)
            A[i3, i2] = sgn * d3(j3, m3, J2, mm, j2, -m2) * d3(j3, n3, J2, nn, j2, -n2)
    return ph * A

def build(beta, betap, jmax2=2):
    B = make_basis(jmax2)
    d = len(B)
    t = {j2: 2 * iv(j2 + 1, beta) / beta for j2 in range(jmax2 + 1)}
    K1 = np.diag([t[b[0]] for b in B])
    Mf = np.zeros((d * d, d * d))
    scale = np.exp(-betap)
    for J2 in range(0, jmax2 + 1):
        if betap == 0.0:
            cJ = 1.0 if J2 == 0 else 0.0
        else:
            cJ = 2 * (J2 + 1) * iv(J2 + 1, betap) / betap
        for M2 in range(-J2, J2 + 1, 2):
            for N2 in range(-J2, J2 + 1, 2):
                A1 = mult_op(B, J2, M2, N2, conj=False)
                A2 = mult_op(B, J2, M2, N2, conj=True)
                Mf += cJ * np.kron(A1, A2)
    Mf *= scale
    K2 = np.kron(K1, K1)
    T = Mf @ K2 @ Mf
    return T, Mf, t, B, d

print("=== sanity 1: at betap = 0, M_f must equal the identity ===")
T, Mf, t, B, d = build(2.5, 0.0, jmax2=2)
print("  dim=%d ; max |M_f - I| = %.3e" % (Mf.shape[0], np.abs(Mf - np.eye(Mf.shape[0])).max()))

print()
print("=== sanity 2: M_f must be symmetric positive ===")
T, Mf, t, B, d = build(2.5, 1.5, jmax2=2)
print("  dim=%d ; max |M_f - M_f^T| = %.3e ; min eig = %.6e"
      % (Mf.shape[0], np.abs(Mf - Mf.T).max(), np.linalg.eigvalsh((Mf + Mf.T) / 2).min()))

print()
print("=== the two-link ratio test ===")
for beta in (2.5,):
    for betap in (0.0, 0.3, 0.8, 1.5, 3.0):
        T, Mf, t, B, d = build(beta, betap, jmax2=2)
        T = (T + T.T) / 2
        ev = np.sort(np.linalg.eigvalsh(T))[::-1]
        rK = float(iv(2, beta) / iv(1, beta))
        print("  beta=%.1f betap=%.1f : top4 = %s"
              % (beta, betap, np.array2string(ev[:4], precision=6)))
        print("      r(T)=%.10f  I2/I1=%.10f  %s"
              % (ev[1] / ev[0], rK, "OK" if ev[1] / ev[0] <= rK + 1e-9 else "FAIL"))
