import numpy as np
from scipy.special import iv

# TELESCOPING:  M[x,y] = sum_{t=|x-y|}^{x+y} (I_{2t} - I_{2t+2}) = I_{2|x-y|} - I_{2(x+y)+2}
# Verify, then write the 2x2 minor in closed form.

def M_tel(x, y, lam):
    return iv(2 * abs(x - y), lam) - iv(2 * (x + y) + 2, lam)

def M_sum(x, y, lam):
    lo, hi = abs(x - y), x + y
    n = int(round(hi - lo))
    return sum(iv(2 * (lo + t), lam) - iv(2 * (lo + t) + 2, lam) for t in range(n + 1))

print("=== verify the telescoping closed form ===")
lam = 1.5
worst = 0.0
for x in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0):
    for y in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0):
        a = M_sum(x, y, lam); b = M_tel(x, y, lam)
        worst = max(worst, abs(a - b) / max(abs(a), 1e-30))
print("  max relative difference (sum vs telescoped) = %.3e" % worst)
print()
print("  =>  (M_f)_{j'j} = I_{2|j-j'|}(lambda) - I_{2(j+j')+2}(lambda)")
print()

# minor in closed form
def minor_closed(j1, j2, k1, k2, lam):
    return M_tel(j1, k1, lam) * M_tel(j2, k2, lam) - M_tel(j1, k2, lam) * M_tel(j2, k1, lam)

print("=== 2x2 minors via the closed form (should all be >= 0 for lambda>0) ===")
js = [k / 2.0 for k in range(0, 9)]
for lam in (-2.0, -0.5, 0.5, 1.5, 3.0):
    mn = np.inf; arg = None; neg = 0; tot = 0
    for i1 in range(len(js)):
        for i2 in range(i1 + 1, len(js)):
            for k1 in range(len(js)):
                for k2 in range(k1 + 1, len(js)):
                    v = minor_closed(js[i1], js[i2], js[k1], js[k2], lam)
                    tot += 1
                    if v < -1e-14: neg += 1
                    if v < mn: mn = v; arg = (js[i1], js[i2], js[k1], js[k2])
    print("  lambda=%+5.1f : total=%4d negative=%3d  min=%+.4e at %s"
          % (lam, tot, neg, mn, arg))

print()
print("=== is the minor factorisable? inspect a few symbolic-looking cases ===")
lam = 1.5
for (j1, j2, k1, k2) in [(0.0, 0.5, 0.0, 0.5), (0.0, 1.0, 0.0, 1.0), (0.0, 1.0, 0.5, 1.0),
                         (0.5, 1.5, 0.5, 1.5), (1.0, 2.0, 1.0, 2.0)]:
    a = M_tel(j1, k1, lam); b = M_tel(j2, k2, lam)
    c = M_tel(j1, k2, lam); d = M_tel(j2, k1, lam)
    print("  (%3.1f,%3.1f;%3.1f,%3.1f): ac=%+.8f bd=%+.8f ad=%+.8f bc=%+.8f  D=%+.6e"
          % (j1, j2, k1, k2, a, b, c, d, a * b - c * d))
