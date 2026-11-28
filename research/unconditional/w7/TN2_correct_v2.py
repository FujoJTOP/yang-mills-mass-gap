import numpy as np
from scipy.special import iv

# CORRECT construction: index the character coefficient by the DOUBLED index m = 2d,
# and sum with step 2 in the doubled index (step 1 in j-units).
#   c2[m] = c_{m/2}
#   M[ia,ib] = sum_{m=|ia-ib|, |ia-ib|+2, ..., ia+ib} c2[m]
# Verified against Haar quadrature in minor_fixed.py to 1e-16.

om = np.linspace(0, np.pi, 40001)
meas = (2 / np.pi) * np.sin(om) ** 2

def ch(j, o):
    s = np.sin(o); num = np.sin((2 * j + 1) * o)
    return np.where(np.abs(s) < 1e-12, (2 * j + 1), num / np.where(np.abs(s) < 1e-12, 1.0, s))

def c2_from_f(fv, J2):
    """character coefficients indexed by m=2j, m=0..J2."""
    return np.array([np.trapezoid(fv * ch(m / 2.0, om) * meas, om) for m in range(J2 + 1)])

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

J2 = 10
b = 1.5
print("=== verify construction against the known closed form c2[m] = e^-b(I_m - I_{m+2}) ===")
c2a = c2_from_f(np.exp(-b * (1 - np.cos(om))), 2 * J2 + 2)
c2b = np.array([np.exp(-b) * (iv(m, b) - iv(m + 2, b)) for m in range(2 * J2 + 3)])
print("  max |quad - closed| = %.3e" % np.abs(c2a - c2b).max())
M = build(c2b, J2)
print("  min 2x2 minor (actual Bessel c) = %+.4e  TN2=%s"
      % (min_minor(M, J2), "YES" if min_minor(M, J2) >= -1e-12 else "NO"))

print()
print("=== conjecture test: positive f  =>  M_f TN2 ? (correct construction) ===")
tests = {
    "f = 1                 ": np.ones_like(om),
    "f = exp(-1.5(1-cos))  ": np.exp(-1.5 * (1 - np.cos(om))),
    "f = exp(+1.5 cos)     ": np.exp(1.5 * np.cos(om)),
    "f = 1 + cos           ": 1 + np.cos(om),
    "f = (1-cos)^2         ": (1 - np.cos(om)) ** 2,
    "f = exp(-5(1-cos)^2)  ": np.exp(-5 * (1 - np.cos(om)) ** 2),
    "f = 1+0.5cos+0.3cos2  ": 1 + 0.5 * np.cos(om) + 0.3 * np.cos(2 * om),
}
for name, fv in tests.items():
    c2 = c2_from_f(fv, 2 * J2 + 2)
    M = build(c2, J2)
    mn = min_minor(M, J2)
    print("  %s : min = %+.4e  TN2=%s" % (name, mn, "YES" if mn >= -1e-12 * max(M.max() ** 2, 1) else "NO"))

print()
print("=== random smooth positive f ===")
rng = np.random.default_rng(17)
fails = 0; tot = 0
for trial in range(15):
    deg = int(rng.integers(1, 4))
    coef = rng.uniform(-2, 2, deg + 1)
    fv = np.exp(sum(coef[k] * np.cos(om) ** k for k in range(deg + 1)))
    c2 = c2_from_f(fv, 2 * J2 + 2)
    M = build(c2, J2)
    mn = min_minor(M, J2)
    tot += 1
    if mn < -1e-10 * max(M.max() ** 2, 1):
        fails += 1
print("  tested %d ; TN2 failures = %d" % (tot, fails))


