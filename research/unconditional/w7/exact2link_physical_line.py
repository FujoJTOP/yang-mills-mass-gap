import numpy as np
from scipy.special import iv
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN

# CRITICAL: in the real gauge model  beta' = beta  (the SAME coupling).
# Along beta'=beta, does r(T) stay bounded away from 1 (=> mass gap), or -> 1 (=> none)?

def d3(j1, m1, j2, m2, j3, m3):
    return float(sN(wigner_3j(Rational(j1, 2), Rational(j2, 2), Rational(j3, 2),
                             Rational(m1, 2), Rational(m2, 2), Rational(m3, 2)), 30))

def make_basis(jmax2):
    B = []
    for j2 in range(0, jmax2 + 1):
        for m2 in range(-j2, j2 + 1, 2):
            for n2 in range(-j2, j2 + 1, 2):
                B.append((j2, m2, n2))
    return B

def mult_op_orth(B, J2, M2, N2, conj=False):
    d = len(B); A = np.zeros((d, d))
    mm, nn = (-M2, -N2) if conj else (M2, N2)
    ph = ((-1.0) ** ((M2 - N2) // 2)) if conj else 1.0
    for i3, (j3, m3, n3) in enumerate(B):
        for i2, (j2, m2, n2) in enumerate(B):
            if m3 != m2 - mm or n3 != n2 - nn: continue
            sgn = (-1.0) ** ((m2 - n2) // 2)
            val = sgn * d3(j3, m3, J2, mm, j2, -m2) * d3(j3, n3, J2, nn, j2, -n2)
            if val == 0.0: continue
            A[i3, i2] = np.sqrt((j3 + 1.0) * (j2 + 1.0)) * val
    return ph * A

JMAX2 = 2
B = make_basis(JMAX2); d = len(B)

def build_Mf(betap):
    Mf = np.zeros((d * d, d * d)); scale = np.exp(-betap)
    for J2 in range(0, JMAX2 + 1):
        cJ = (1.0 if J2 == 0 else 2 * (J2 + 1) * iv(J2 + 1, betap) / betap) if betap > 0 else (1.0 if J2 == 0 else 0.0)
        for M2 in range(-J2, J2 + 1, 2):
            for N2 in range(-J2, J2 + 1, 2):
                A1 = mult_op_orth(B, J2, M2, N2, False)
                A2 = mult_op_orth(B, J2, M2, N2, True)
                Mf += cJ * np.kron(A1, A2)
    return Mf * scale

print("=== along the physical line beta' = beta ===")
print("%8s %14s %16s %14s" % ("beta", "r(T)", "lambda1", "beta*(1-r)"))
for beta in (0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0):
    Mf = build_Mf(beta)
    tj = {j2: 2 * iv(j2 + 1, beta) / beta for j2 in range(JMAX2 + 1)}
    K = np.diag([tj[b1[0]] * tj[b2[0]] for b1 in B for b2 in B])
    T = Mf @ K @ Mf
    T = (T + T.T) / 2
    ev = np.sort(np.linalg.eigvalsh(T))[::-1]
    r = ev[1] / ev[0]
    print("%8.1f %14.10f %16.8e %14.6f" % (beta, r, ev[0], beta * (1 - r)))

print()
print("=== control: same but with V removed (betap = 0) ===")
print("%8s %14s" % ("beta", "r(T) at betap=0 = I2/I1"))
for beta in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
    print("%8.1f %14.10f" % (beta, float(iv(2, beta) / iv(1, beta))))
