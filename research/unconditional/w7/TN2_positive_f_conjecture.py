import numpy as np

# CONJECTURE:  f >= 0 on SU(2)  ==>  M_f (multiplication by f, in the character basis)
#              is totally nonnegative (TN2).
# Test: sample random positive f (as a function of omega), compute c_d = int f chi_d dU,
#       build M, check all 2x2 minors >= 0.

om = np.linspace(0, np.pi, 20001)
meas = (2 / np.pi) * np.sin(om) ** 2

def ch(j, o):
    s = np.sin(o); num = np.sin((2 * j + 1) * o)
    return np.where(np.abs(s) < 1e-12, (2 * j + 1), num / np.where(np.abs(s) < 1e-12, 1.0, s))

def chars(fvals, J):
    return np.array([np.trapezoid(fvals * ch(k / 2.0, om) * meas, om) for k in range(J + 1)])

def build_M(c, J):
    M = np.zeros((J + 1, J + 1))
    for ia in range(J + 1):
        for ib in range(J + 1):
            lo = int(round(abs(ia - ib) / 2.0))
            hi = int(round((ia + ib) / 2.0))
            M[ia, ib] = sum(c[lo:hi + 1])
    return M

def min_minor(M, J):
    mn = np.inf
    for i in range(J):
        for j in range(i + 1, J + 1):
            for k in range(J):
                for l in range(k + 1, J + 1):
                    mn = min(mn, M[i, k] * M[j, l] - M[i, l] * M[j, k])
    return mn

rng = np.random.default_rng(17)
J = 8

print("=== test: positive f  =>  M_f TN2 ? ===")
tests = {
    "f = 1                    ": np.ones_like(om),
    "f = exp(-1.5(1-cos))     ": np.exp(-1.5 * (1 - np.cos(om))),
    "f = exp(+1.5 cos)        ": np.exp(1.5 * np.cos(om)),
    "f = 1 + cos              ": 1 + np.cos(om),
    "f = (1-cos)^2            ": (1 - np.cos(om)) ** 2,
    "f = exp(-5(1-cos)^2)     ": np.exp(-5 * (1 - np.cos(om)) ** 2),
    "f = 1+0.5cos+0.3cos2     ": 1 + 0.5 * np.cos(om) + 0.3 * np.cos(2 * om),
}
for name, fv in tests.items():
    if fv.min() < -1e-14:
        print("  %s : f NOT >= 0, skip" % name); continue
    c = chars(fv, J)
    M = build_M(c, J)
    mn = min_minor(M, J)
    print("  %s : min 2x2 minor = %+.4e  TN2=%s" % (name, mn, "YES" if mn >= -1e-12 * M.max() ** 2 else "NO"))

print()
print("=== random positive f (smooth, as polynomials in cos) ===")
fails = 0; tot = 0
for trial in range(20):
    deg = rng.integers(1, 4)
    coef = rng.uniform(-2, 2, deg + 1)
    fv = sum(coef[k] * np.cos(om) ** k for k in range(deg + 1))
    fv = np.exp(fv)                      # ensure positive
    c = chars(fv, J)
    M = build_M(c, J)
    mn = min_minor(M, J)
    tot += 1
    if mn < -1e-10 * M.max() ** 2:
        fails += 1
print("  tested %d random positive f ; TN2 failures = %d" % (tot, fails))

print()
print("=== control: f NOT >= 0 (sign-changing) ===")
fails = 0; tot = 0
for trial in range(20):
    deg = rng.integers(1, 4)
    coef = rng.uniform(-2, 2, deg + 1)
    fv = sum(coef[k] * np.cos(om) ** k for k in range(deg + 1))
    if fv.min() >= -1e-14:
        continue
    c = chars(fv, J)
    M = build_M(c, J)
    mn = min_minor(M, J)
    tot += 1
    if mn < -1e-10 * M.max() ** 2:
        fails += 1
print("  tested %d sign-changing f ; TN2 failures = %d" % (tot, fails))
