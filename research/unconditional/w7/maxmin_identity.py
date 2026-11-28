import numpy as np
from scipy.special import iv

# Derived identity to verify:
#   (M_f)_{j'j} = (2/pi) * int_0^pi e^{lambda cos th} * sin((2*max(j,j')+1)th) * sin((2*min(j,j')+1)th) dth
# Derivation used:  sum_{t=p}^{q} sin((2t+1)th) = sin((p+q+1)th) sin((q-p+1)th)/sin(th)
#   with p=|j-j'|, q=j+j'  =>  p+q = 2*max, q-p = 2*min.

def M_bessel(jp, j, lam):
    a = int(round(2 * jp)); b = int(round(2 * j))
    lo, hi = abs(a - b), a + b
    return sum(np.exp(-0.0) * (iv(m, lam) - iv(m + 2, lam)) for m in range(lo, hi + 1, 2))

def M_maxmin(jp, j, lam, nq=200001):
    th = np.linspace(0, np.pi, nq)
    mx, mn = max(jp, j), min(jp, j)
    return (2 / np.pi) * np.trapezoid(np.exp(lam * np.cos(th)) *
                                      np.sin((2 * mx + 1) * th) * np.sin((2 * mn + 1) * th), th)

print("=== verify the max-min identity ===")
lam = 1.5
worst = 0.0
for jp in (0.0, 0.5, 1.0, 1.5, 2.0):
    for j in (0.0, 0.5, 1.0, 1.5, 2.0):
        a = M_bessel(jp, j, lam); b = M_maxmin(jp, j, lam)
        worst = max(worst, abs(a - b) / max(abs(a), 1e-30))
print("  lambda=%.1f : max relative difference (Bessel vs max-min) = %.3e" % (lam, worst))
print()

print("=== why the lambda-threshold: where does e^{lambda cos} concentrate? ===")
print("  lambda>0 : weight peaks at theta=0,  where sin((2j+1)th) ~ (2j+1)th  (no sign alternation)")
print("  lambda<0 : weight peaks at theta=pi, where sin((2j+1)th) ~ (2j+1)(pi-th)*(-1)^j (alternating!)")
print()
print("=== the 2x2 minor as a double integral over the antisymmetric determinants ===")
print("  Delta = (2/pi)^2 int int e^{lam(cos th + cos ph)} D_ab(th,ph) D_cd(th,ph) dth dph")
print("  with D_ab = u_a(th)u_b(ph) - u_a(ph)u_b(th),  u_j = sin((2j+1)*.)")
print()
print("  Check: is D_ab(th,ph) of fixed sign on the triangle th+ph<pi vs >pi ?")
for (a, b) in [(0, 1), (0, 2), (1, 2)]:
    th = np.linspace(0.01, np.pi - 0.01, 400)
    TH, PH = np.meshgrid(th, th, indexing='ij')
    Ua = lambda x: np.sin((2 * a + 1) * x); Ub = lambda x: np.sin((2 * b + 1) * x)
    D = Ua(TH) * Ub(PH) - Ua(PH) * Ub(TH)
    low = D[TH + PH < np.pi - 1e-3]; high = D[TH + PH > np.pi + 1e-3]
    print("    (a,b)=(%d,%d): D range on th+ph<pi : [%+.4f, %+.4f] ; on th+ph>pi : [%+.4f, %+.4f]"
          % (a, b, low.min(), low.max(), high.min(), high.max()))
