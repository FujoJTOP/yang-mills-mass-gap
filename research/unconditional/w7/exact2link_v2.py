import numpy as np
from scipy.special import iv

# FIXED exact SU(2) two-link computation (Peter-Weyl truncated, no symmetry assumptions).
#   basis per link: |j,m,n>, 2j in {0,1,2,3} (j <= 3/2), dim = 1+4+9+16 = 30
#   K1 diagonal: t_j = 2*I_{2j+1}(beta)/beta
#   f(U1 U2^-1) = exp(-betap*(1 - cos(theta/2))),  expanded as f = sum_j c_j chi_j
#       c_j = int f chi_j dU   (chi_j orthonormal, Haar dU = (2/pi) sin^2(theta/2) dtheta)
#   chi_j(U1 U2^-1) = sum_{m,n} D^j_{mn}(U1) conj(D^j_{mn}(U2))
#   => M_f = sum_j c_j sum_{m,n} M1[D^j_{mn}] (x) M2[conj D^j_{mn}]
#   multiplication matrix elements from 3j symbols:
#       int D^{j1}_{m1n1} D^{j2}_{m2n2} conj(D^{j3}_{m3n3}) dU
#         = 3j(j1 j2 j3; m1 m2 -m3) * 3j(j1 j2 j3; n1 n2 -n3)     (up to a phase from conj)
#   NOTE: conj(D^j_{mn}) = (-1)^(m-n) D^j_{-m,-n}

def charval(j2, theta):
    # j2 = 2j ; chi_j = sin((2j+1) theta/2)/sin(theta/2)
    x = theta / 2.0
    s = np.sin(x)
    num = np.sin((j2 + 1.0) * x)
    out = np.where(np.abs(s) < 1e-12, (j2 + 1.0), num / np.where(np.abs(s) < 1e-12, 1.0, s))
    return out

def ccoeffs(betap, jmax2, nq=2001):
    th = np.linspace(0, np.pi, nq)
    q = np.cos(th / 2.0)
    f = np.exp(-betap * (1.0 - q))
    meas = (2.0 / np.pi) * np.sin(th / 2.0) ** 2
    out = []
    for j2 in range(0, jmax2 + 1):
        ch = charval(j2, th)
        out.append(np.trapezoid(f * ch * meas, th))
    return np.array(out)

print("=== character coefficients of f = exp(-betap(1-cos))  (must be POSITIVE) ===")
for betap in (0.5, 1.5, 3.0):
    c = ccoeffs(betap, 12)
    print("  betap=%.1f : c_j (j2=0..6) = %s" % (betap, np.array2string(c[:7], precision=6)))
print()
print("  cross-check: f(U) = sum_j c_j chi_j(U) must reproduce f at theta")
for betap in (1.5,):
    c = ccoeffs(betap, 40)
    th = np.array([0.0, 0.7, 1.4, 2.1, 2.8])
    f_true = np.exp(-betap * (1.0 - np.cos(th / 2.0)))
    f_ser = sum(c[j2] * charval(j2, th) for j2 in range(len(c)))
    print("   theta :", th)
    print("   f_true:", np.array2string(f_true, precision=6))
    print("   f_ser :", np.array2string(f_ser, precision=6))

