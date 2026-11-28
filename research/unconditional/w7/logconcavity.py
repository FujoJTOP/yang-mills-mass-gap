import numpy as np
from scipy.special import iv

# Is the actual c_d = e^{-b}(I_{2d}(b) - I_{2d+2}(b)) log-concave?
# (log r_d = c_d^2 / (c_{d-1} c_{d+1}) >= 1)
# And: does log-concavity of c_d  <==>  M_f TN2 ?

def cd(d, b):
    return np.exp(-b) * (iv(2 * d, b) - iv(2 * d + 2, b))

print("=== log-concavity ratio r_d = c_d^2/(c_{d-1}c_{d+1})  (>=1 means log-concave) ===")
for b in (0.5, 1.5, 3.0, 8.0):
    c = np.array([cd(k, b) for k in range(25)])
    r = np.array([c[k] ** 2 / (c[k - 1] * c[k + 1]) for k in range(1, len(c) - 1)])
    print("  b=%.1f : min r_d = %.6f  (d=1..%d)   log-concave=%s"
          % (b, r.min(), len(c) - 2, "YES" if r.min() >= 1 - 1e-12 else "NO"))
    print("        r_d[:6] = %s" % np.array2string(r[:6], precision=6))

print()
print("=== is c_d also PF2 (Toeplitz TN2) : c_{d}c_{d+2} - c_{d+1}^2 ? (sign) ===")
for b in (0.5, 1.5, 3.0):
    c = np.array([cd(k, b) for k in range(25)])
    pf = np.array([c[k] * c[k + 2] - c[k + 1] ** 2 for k in range(len(c) - 2)])
    print("  b=%.1f : c_d c_{d+2} - c_{d+1}^2 : min=%+.3e max=%+.3e  => sign %s"
          % (b, pf.min(), pf.max(), "always <=0 (log-concave)" if pf.max() <= 1e-18 else "MIXED"))

print()
print("=== log-concavity of the Bessel ratios I_{n+1}/I_n  (known: decreasing in n) ===")
for b in (1.5,):
    I = np.array([iv(n, b) for n in range(0, 20)])
    ratio = I[1:] / I[:-1]
    print("  b=%.1f : I_{n+1}/I_n decreasing = %s" % (b, np.all(np.diff(ratio) <= 1e-15)))
    print("        values[:8] = %s" % np.array2string(ratio[:8], precision=6))
