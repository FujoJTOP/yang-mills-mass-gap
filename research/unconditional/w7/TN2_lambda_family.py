import numpy as np
from scipy.special import iv

# Sharper conjecture: M_f is TN2 for f = e^{lambda q} (q = cos omega), for ALL lambda.
#   For f = e^{lambda q} the character coefficients are  c_m = 2(m+1) I_{m+1}(lambda)/lambda
#   (m = 2d).  The constant factor of e^{-b(1-q)} does not affect TN2.
# Also map which simple families work / fail.

om = np.linspace(0, np.pi, 40001)
meas = (2 / np.pi) * np.sin(om) ** 2
q = np.cos(om)

def ch(j, o):
    s = np.sin(o); num = np.sin((2 * j + 1) * o)
    return np.where(np.abs(s) < 1e-12, (2 * j + 1), num / np.where(np.abs(s) < 1e-12, 1.0, s))

def c2_from_f(fv, M2):
    return np.array([np.trapezoid(fv * ch(m / 2.0, om) * meas, om) for m in range(M2 + 1)])

def build(c2, J2):
    M = np.zeros((J2 + 1, J2 + 1))
    for ia in range(J2 + 1):
        for ib in range(J2 + 1):
            lo, hi = abs(ia - ib), ia + ib
            M[ia, ib] = sum(c2[m] for m in range(lo, hi + 1, 2) if m < len(c2))
    return M

def min_minor(M, J2):
    mn = np.inf
    for i in range(J2):
        for j in range(i + 1, J2 + 1):
            for k in range(J2):
                for l in range(k + 1, J2 + 1):
                    mn = min(mn, M[i, k] * M[j, l] - M[i, l] * M[j, k])
    return mn

J2 = 8
M2 = 2 * J2 + 2
print("=== f = e^{lambda q} : TN2 for which lambda? ===")
for lam in (-6, -4, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 6):
    fv = np.exp(lam * q)
    c2 = c2_from_f(fv, M2)
    M = build(c2, J2)
    mn = min_minor(M, J2)
    print("  lambda=%+5.1f : min minor = %+.4e  TN2=%s"
          % (lam, mn, "YES" if mn >= -1e-11 * max(M.max() ** 2, 1) else "NO"))

print()
print("=== f = q^k (monomials) ===")
for k in (0, 1, 2, 3):
    fv = q ** k
    if fv.min() < 0: fv = fv + 1e-9
    c2 = c2_from_f(fv, M2)
    M = build(c2, J2)
    mn = min_minor(M, J2)
    print("  f=q^%d : min minor = %+.4e  TN2=%s" % (k, mn, "YES" if mn >= -1e-11 * max(M.max() ** 2, 1) else "NO"))

print()
print("=== f = (1+aq)^k  (generalised) ===")
for k in (1, 2, 3):
    for a in (0.3, 0.6):
        fv = (1 + a * q) ** k
        c2 = c2_from_f(fv, M2)
        M = build(c2, J2)
        mn = min_minor(M, J2)
        print("  k=%d a=%.1f : min = %+.4e  TN2=%s" % (k, a, mn, "YES" if mn >= -1e-11 * max(M.max() ** 2, 1) else "NO"))
