import numpy as np
from scipy.special import iv

# DECISIVE CHECK: is M(lambda) = M(-lambda) ?
# Analytic: (M_f)_{j'j} = I_{2|j-j'|}(lambda) - I_{2(j+j')+2}(lambda), all orders EVEN
#           and I_n(-z) = (-1)^n I_n(z)  =>  M(-lambda) = M(lambda).
# Also: M[x,y] = (2/pi) int_0^pi e^{lambda cos th} u_x u_y dth, and theta->pi-theta
#       maps cos -> -cos while u_j(pi-th) = u_j(th). So M(lambda)=M(-lambda) again.

print("=== scipy.special.iv for negative argument, integer order ===")
for n in (0, 1, 2, 3, 4):
    print("  iv(%d, -1.5) = %+.10f   (-1)^n iv(%d, 1.5) = %+.10f"
          % (n, iv(n, -1.5), n, (-1) ** n * iv(n, 1.5)))

print()
print("=== M_tel(lambda) vs M_tel(-lambda) ===")
def M_tel(x, y, lam):
    return iv(2 * abs(x - y), lam) - iv(2 * (x + y) + 2, lam)

worst = 0.0
for x in (0.0, 0.5, 1.0, 1.5, 2.0):
    for y in (0.0, 0.5, 1.0, 1.5, 2.0):
        a = M_tel(x, y, 1.5); b = M_tel(x, y, -1.5)
        worst = max(worst, abs(a - b))
print("  max |M(1.5) - M(-1.5)| = %.3e" % worst)

print()
print("=== quadrature version: int e^{lam cos} u_x u_y ===")
om = np.linspace(0, np.pi, 200001)
meas = (2 / np.pi) * np.sin(om) ** 2
def ch(j, o):
    s = np.sin(o); num = np.sin((2 * j + 1) * o)
    return np.where(np.abs(s) < 1e-12, (2 * j + 1), num / np.where(np.abs(s) < 1e-12, 1.0, s))
worst = 0.0
for lam in (0.5, 1.5, 3.0):
    for (x, y) in [(0.0, 0.0), (0.5, 0.0), (1.0, 0.5), (1.5, 0.5)]:
        a = np.trapezoid(np.exp(lam * np.cos(om)) * ch(x, om) * ch(y, om) * meas, om)
        b = np.trapezoid(np.exp(-lam * np.cos(om)) * ch(x, om) * ch(y, om) * meas, om)
        worst = max(worst, abs(a - b))
print("  max |M(+lam) - M(-lam)| over lam in {0.5,1.5,3.0} = %.3e" % worst)

print()
print("=== so the lambda-SIGN threshold must be an artifact. Recheck TN2 at +-lam ===")
js = [k / 2.0 for k in range(0, 9)]
def minor(j1, j2, k1, k2, lam):
    return M_tel(j1, k1, lam) * M_tel(j2, k2, lam) - M_tel(j1, k2, lam) * M_tel(j2, k1, lam)
for lam in (-3.0, -1.5, -0.5, 0.5, 1.5, 3.0):
    mn = np.inf; neg = 0; tot = 0
    for i1 in range(len(js)):
        for i2 in range(i1 + 1, len(js)):
            for k1 in range(len(js)):
                for k2 in range(k1 + 1, len(js)):
                    v = minor(js[i1], js[i2], js[k1], js[k2], lam)
                    tot += 1
                    if v < -1e-12: neg += 1
                    if v < mn: mn = v
    print("  lam=%+5.1f : neg=%4d/%d  min=%+.4e" % (lam, neg, tot, mn))
