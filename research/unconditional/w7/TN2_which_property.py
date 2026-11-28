import numpy as np
from scipy.special import iv

# Which property of c_d makes M_f totally nonnegative (TN2)?
#   (M_f)_{j'j} = sum_{d=|j-j'|}^{j+j'} c_d
# Test M_f TN2 for several families of c_d:
#   (a) c_d = e^-b (I_2d - I_2d+2)      [the actual coefficients]  (positive, decreasing)
#   (b) c_d = r^d, 0<r<1                (geometric, log-concave)
#   (c) c_d = 1/(d+1)^p                 (algebraic)
#   (d) c_d = 1 for d<=D else 0         (step)
#   (e) random positive decreasing
#   (f) random positive NON-monotone

def build(c, J):
    M = np.zeros((J + 1, J + 1))
    for a in range(J + 1):
        for b in range(J + 1):
            ja, jb = a / 2.0, b / 2.0
            lo, hi = abs(ja - jb), ja + jb
            s = 0.0
            d = lo
            while d <= hi + 1e-12:
                k = int(round(d))
                if 0 <= k < len(c):
                    s += c[k]
                d += 1.0
            M[a, b] = s
    return M

def min_minor(M, J):
    mn = np.inf
    for i in range(J):
        for j in range(i + 1, J + 1):
            for k in range(J):
                for l in range(k + 1, J + 1):
                    mn = min(mn, M[i, k] * M[j, l] - M[i, l] * M[j, k])
    return mn

J = 10
b = 1.5
fams = {}
fams["(a) e^-b(I_2d-I_2d+2)"] = np.array([np.exp(-b) * (iv(2 * d, b) - iv(2 * d + 2, b)) for d in range(J + 1)])
fams["(b) 0.5^d            "] = 0.5 ** np.arange(J + 1)
fams["(c) 1/(d+1)^3        "] = 1.0 / (np.arange(J + 1) + 1.0) ** 3
fams["(d) step D=4         "] = np.array([1.0 if d <= 4 else 0.0 for d in range(J + 1)])
rng = np.random.default_rng(5)
v = np.sort(rng.uniform(0.01, 1.0, J + 1))[::-1]
fams["(e) random decreasing"] = v
w = rng.uniform(0.01, 1.0, J + 1)
fams["(f) random non-monot "] = w

print("=== M_f TN2 for different c_d families (J=10) ===")
for name, c in fams.items():
    M = build(c, J)
    mn = min_minor(M, J)
    mono = np.all(np.diff(c) <= 1e-15)
    print("  %s : min 2x2 minor = %+.3e   decreasing(c)=%s   TN2=%s"
          % (name, mn, mono, "YES" if mn >= -1e-10 else "NO"))

print()
print("=== is the indicator family I_d itself TN2? ===")
for d in (0, 1, 2, 3):
    c = np.array([1.0 if k == d else 0.0 for k in range(J + 1)])
    M = build(c, J)
    mn = min_minor(M, J)
    print("  I_%d : min 2x2 minor = %+.3e  TN2=%s" % (d, mn, "YES" if mn >= -1e-10 else "NO"))
