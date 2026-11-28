import numpy as np
from scipy.special import iv
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
import itertools, functools

# EXACT two-link SU(2) computation in the Peter-Weyl basis |j,m,n>, j <= jmax.
#   f(u1,u2) = exp(-betap*(1 - <q1,q2>)),  <q1,q2> = (1/2) tr(u1 u2^-1) = cos(omega)
#   f = exp(-betap) * sum_J c_J chi_J ,   c_J = 2(2J+1) I_{2J+1}(betap)/betap
#   chi_J(u1 u2^-1) = sum_{M,N} D^J_{MN}(u1) conj(D^J_{MN}(u2))
#   => M_f = exp(-betap) sum_J c_J sum_{M,N} M1[D^J_{MN}] (x) M2[conj D^J_{MN}]
#   matrix elements from  int D^{j1}_{m1n1} D^{j2}_{m2n2} D^{j3}_{m3n3} dU
#                                  = 3j(j1 j2 j3; m1 m2 m3) * 3j(j1 j2 j3; n1 n2 n3)
#   T = M_f K M_f,  K = K1 (x) K1,  K1 diagonal with t_j = 2 I_{2j+1}(beta)/beta

def d3(j1, m1, j2, m2, j3, m3):
    return float(sN(wigner_3j(Rational(j1), Rational(m1), Rational(j2), Rational(m2),
                             Rational(j3), Rational(m3)), 25))

def idx(j2, m2, n2):
    return (j2, m2, n2)

def make_basis(jmax2):
    """jmax2 = 2*jmax. Returns list of (2j, 2m, 2n)."""
    B = []
    for j2 in range(0, jmax2 + 1):
        for m2 in range(-j2, j2 + 1, 2):
            for n2 in range(-j2, j2 + 1, 2):
                B.append((j2, m2, n2))
    return B

def mult_first(B, J2, M2, N2):
    """M1[D^J_{MN}] acting on the first factor, in basis B (2j,2m,2n ints)."""
    d = len(B)
    A = np.zeros((d, d))
    # element from |j2,m2,n2> to |j3,m3,n3>:
    # int D^{j3}_{m3 n3} D^J_{M N} conj(D^{j2}_{m2 n2}) dU
    # conj(D^{j2}_{m2 n2}) = (-1)^{m2-n2} D^{j2}_{-m2,-n2}
    for i3, (j3, m3, n3) in enumerate(B):
        for i2, (j2, m2, n2) in enumerate(B):
            if m3 != M2 + (-m2) or n3 != N2 + (-n2):
                continue
            sgn = (-1) ** ((m2 - n2) // 2)
            a = d3(j3, m3, J2, M2, j2, -m2)
            b = d3(j3, n3, J2, N2, j2, -n2)
            A[i3, i2] = sgn * a * b
    return A

def build(beta, betap, jmax2=2):
    B = make_basis(jmax2)
    d = len(B)
    t = {j2: 2 * iv(j2 + 1, beta) / beta for j2 in range(jmax2 + 1)}
    # K = K1 (x) K1 diagonal
    K1 = np.diag([t[b[0]] for b in B])
    # M_f
    Mf = np.zeros((d * d, d * d))
    scale = np.exp(-betap)
    for J2 in range(0, jmax2 + 1):
        cJ = 2 * (J2 + 1) * iv(J2 + 1, betap) / betap
        for M2 in range(-J2, J2 + 1, 2):
            for N2 in range(-J2, J2 + 1, 2):
                A = mult_first(B, J2, M2, N2)
                As = A * ((-1) ** ((M2 - N2) // 2))     # conj factor
                Blk = np.kron(As, A)
                Mf += cJ * Blk
    Mf *= scale
    K2 = np.kron(K1, K1)
    T = Mf @ K2 @ Mf
    return T, t, B

print("building exact 2-link SU(2) model, j <= 1 (dim 14 per link, 196 total) ...")
for beta, betap in [(2.5, 1.5)]:
    T, t, B = build(beta, betap, jmax2=2)
    T = (T + T.T) / 2
    ev = np.sort(np.linalg.eigvalsh(T))[::-1]
    rK = float(iv(2, beta) / iv(1, beta))
    print("  beta=%.2f betap=%.2f : dim=%d" % (beta, betap, T.shape[0]))
    print("    top5 eigenvalues = %s" % np.array2string(ev[:5], precision=8))
    print("    r(T) = t2/t1 = %.10f    I2/I1 = %.10f    %s"
          % (ev[1] / ev[0], rK, "OK" if ev[1] / ev[0] <= rK + 1e-10 else "FAIL"))
