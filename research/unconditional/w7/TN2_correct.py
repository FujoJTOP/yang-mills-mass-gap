import numpy as np
from scipy.special import iv

# CORRECT construction (round 13 was right; round 14 had a banker's-rounding bug):
#   (M_f)_{j'j} = sum_{d=|j-j'|}^{j+j'} c_d   with d in j-units (step 1)
#   c_d = e^{-b}(I_{2d}(b) - I_{2d+2}(b))
# Index the matrix by a = 2j' (integer), so d (j-units) takes values |a-b|/2 .. (a+b)/2 step 1.

def build(cfun, Jf, b):
    """cfun(d, b) -> c_d,  d in j-units (may be half-integer)."""
    M = np.zeros((Jf + 1, Jf + 1))
    for ia in range(Jf + 1):
        for ib in range(Jf + 1):
            ja, jb = ia / 2.0, ib / 2.0
            lo, hi = abs(ja - jb), ja + jb
            n = int(round((hi - lo)))            # number of steps, all integers in j-units
            s = 0.0
            for t in range(n + 1):
                s += cfun(lo + t, b)
            M[ia, ib] = s
    return M

def min_minor(M, Jf):
    mn = np.inf; arg = None
    for i in range(Jf):
        for j in range(i + 1, Jf + 1):
            for k in range(Jf):
                for l in range(k + 1, Jf + 1):
                    dt = M[i, k] * M[j, l] - M[i, l] * M[j, k]
                    if dt < mn: mn = dt; arg = (i, j, k, l)
    return mn, arg

b = 1.5; Jf = 10

print("=== sanity: entries vs exact Haar integral ===")
om = np.linspace(0, np.pi, 400001); meas = (2 / np.pi) * np.sin(om) ** 2
f = np.exp(-b * (1 - np.cos(om)))
def charval(j, o):
    s = np.sin(o); num = np.sin((2 * j + 1) * o)
    return np.where(np.abs(s) < 1e-12, (2 * j + 1), num / np.where(np.abs(s) < 1e-12, 1.0, s))
M = build(lambda d, bb: np.exp(-bb) * (iv(2 * d, bb) - iv(2 * d + 2, bb)), Jf, b)
for (jp, j) in [(0.0, 0.0), (0.5, 0.0), (0.5, 0.5), (1.0, 0.5), (1.5, 1.0)]:
    q = np.trapezoid(charval(jp, om) * f * charval(j, om) * meas, om)
    a1, a2 = int(round(2 * jp)), int(round(2 * j))
    print("  (j'=%3.1f,j=%3.1f) quad=%.10f  build=%.10f  diff=%.2e" % (jp, j, q, M[a1, a2], abs(q - M[a1, a2])))

print()
print("=== TN2 for different c_d families (correct construction) ===")
fams = {
    "(a) e^-b(I_2d-I_2d+2) ": lambda d, bb: np.exp(-bb) * (iv(2 * d, bb) - iv(2 * d + 2, bb)),
    "(b) 0.5^d             ": lambda d, bb: 0.5 ** d,
    "(c) 1/(d+1)^3         ": lambda d, bb: 1.0 / (d + 1.0) ** 3,
    "(d) step D=4          ": lambda d, bb: 1.0 if d <= 4 else 0.0,
}
rng = np.random.default_rng(5)
vl = np.sort(rng.uniform(0.01, 1.0, 30))[::-1]
fams["(e) random decreasing"] = lambda d, bb: vl[int(round(d))]
wl = rng.uniform(0.01, 1.0, 30)
fams["(f) random non-monot "] = lambda d, bb: wl[int(round(d))]

for name, cf in fams.items():
    Mn = build(cf, Jf, b)
    mn, arg = min_minor(Mn, Jf)
    print("  %s : min 2x2 minor = %+.4e at %s  TN2=%s" % (name, mn, arg, "YES" if mn >= -1e-10 else "NO"))
