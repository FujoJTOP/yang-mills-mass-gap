import numpy as np
from scipy.special import iv
from sympy.physics.wigner import wigner_3j
from sympy import Rational, N as sN
import itertools, functools

# EXACT (Peter-Weyl truncated) two-link SU(2) computation.
# Basis per link: |j,m,n>,  j in {0,1/2,1,...,jmax},  m,n in -j..j.
#   dim = sum_j (2j+1)^2
# Single-link kernel (convolution): K1 diagonal, eigenvalue t_j = 2*I_{2j+1}(beta)/beta
# "spatial" weight f = exp(-betap*(1 - <q1,q2>)),  <q1,q2> = (1/2) tr(U1 U2^dagger)
#   expand f as a function of the relative element:  f = sum_j c_j chi_j
#   with c_j = (2j+1) * int f(U) conj(chi_j(U)) dU   (1D quadrature in theta)
#   then M_f acts on L^2(G x G) via  chi_j(U1 U2^-1) = sum_{m,n} D^j_{mn}(U1) conj(D^j_{mn}(U2))
#   => M_f = sum_j c_j sum_{m,n} M[D^j_{mn} on u1] (x) M[conj(D^j_{mn}) on u2]
#   and multiplication matrix elements come from 3j symbols.

def basis(jmax):
    js = [Rational(2 * k, 2) for k in range(0, 2 * jmax + 1)]
    idx = []
    for j in js:
        ji = int(2 * j)
        for m in range(-ji, ji + 1, 2):
            for n in range(-ji, ji + 1, 2):
                idx.append((ji, m, n))     # store 2j, 2m, 2n as ints
    return idx

def threej(j1, j2, j3, m1, m2, m3):
    return float(sN(wigner_3j(Rational(j1, 2), Rational(j2, 2), Rational(j3, 2),
                             Rational(m1, 2), Rational(m2, 2), Rational(m3, 2)), 20))

def integ(j1, m1, n1, j2, m2, n2, j3, m3, n3):
    # int D^{j1}_{m1n1} D^{j2}_{m2n2} D^{j3}_{m3n3} dU
    # = (2j3+1)^{-1} * 3j(j1 j2 j3; m1 m2 m3) * 3j(j1 j2 j3; n1 n2 n3)
    if m1 + m2 + m3 != 0 or n1 + n2 + n3 != 0:
        return 0.0
    a = threej(j1, j2, j3, m1, m2, m3)
    b = threej(j1, j2, j3, n1, n2, n3)
    return (2.0 * (j3 / 2.0) + 1.0) ** (-1) * a * b

print("building basis ...")
jmax = 3            # j <= 3/2  => dim 1+4+9+16 = 30
B = basis(jmax)
d = len(B)
print("  per-link dim =", d, " ; two-link dim =", d * d)

# characters chi_j(U) at theta: chi_j = sin((2j+1)theta/2)/sin(theta/2), with U angle theta
def fj_theta(theta, j2, beta):
    # f(U) = exp(-betap*(1 - cos(theta/2)))
    return None

def charval(j2, theta):
    x = theta / 2.0
    num = np.sin((j2 + 1) * x)          # chi_j = sin((2j+1)theta/2)/sin(theta/2); (2j+1)/2 = j2+1 ... careful
    den = np.sin(x)
    return np.where(np.abs(den) < 1e-14, (j2 + 1.0), num / np.where(np.abs(den) < 1e-14, 1, den))

def ccoeffs(betap, jmax, nq=400):
    # c_j = (2j+1) * int f chi_j dU,  dU = (2/pi) sin^2(theta/2) dtheta on [0,pi]
    th = np.linspace(0, np.pi, nq)
    w = np.full_like(th, (np.pi / (nq - 1)))
    q = np.cos(th / 2.0)                       # (1/2)tr U
    f = np.exp(-betap * (1.0 - q))
    meas = (2.0 / np.pi) * np.sin(th / 2.0) ** 2
    out = []
    for j2 in range(0, 2 * jmax + 1):
        ch = charval(j2, th)
        val = np.sum(f * ch * meas * w) * (j2 + 1.0)     # (2j+1) = j2+1
        out.append(val)
    return np.array(out)

print("computing character coefficients of f ...")
for betap in (0.5, 1.5, 3.0):
    c = ccoeffs(betap, jmax)
    print("  betap=%.1f : c_j = %s" % (betap, np.array2string(c, precision=6)))
