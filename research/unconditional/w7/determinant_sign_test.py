import numpy as np

# Crux test for the proof of TN2 via the max-min identity:
#   M[x,y] = (2/pi) int w(th) u_max u_min dth,  u_j = sin((2j+1)*),  w = e^{lam cos} > 0
#   =>  Delta = (2/pi)^2 int int w(th)w(ph) D_ab(th,ph) D_cd(th,ph) dth dph
#      with D_ab = u_a(th)u_b(ph) - u_a(ph)u_b(th)
# If D_ab * D_cd >= 0 POINTWISE for a<b, c<d, then Delta >= 0 immediately.

def u(j, x):
    return np.sin((2 * j + 1) * x)

th = np.linspace(1e-6, np.pi - 1e-6, 601)
TH, PH = np.meshgrid(th, th, indexing='ij')

def D(a, b):
    return u(a, TH) * u(b, PH) - u(a, PH) * u(b, TH)

print("=== pointwise sign of D_ab * D_cd  for a<b, c<d ===")
pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3), (0, 5), (2, 5)]
fails = 0; tot = 0
worst = np.inf
for (a, b) in pairs:
    for (c, d) in pairs:
        prod = D(a, b) * D(c, d)
        mn = prod.min()
        tot += 1
        if mn < -1e-9:
            fails += 1
        worst = min(worst, mn)
print("  pairs tested = %d ; violations (prod < 0 somewhere) = %d ; global min = %+.4e"
      % (tot, fails, worst))

print()
print("=== detail for a few ===")
for (a, b, c, d) in [(0, 1, 0, 1), (0, 1, 1, 2), (0, 1, 0, 2), (0, 2, 1, 3), (0, 1, 2, 3)]:
    prod = D(a, b) * D(c, d)
    print("  (%d,%d)x(%d,%d): min=%+.4e max=%+.4e" % (a, b, c, d, prod.min(), prod.max()))

print()
print("=== if the pointwise test fails, check the SIGNED integral instead ===")
# weighted double integral with w(th)w(ph), w = e^{lam cos}
for lam in (0.5, 1.5, 3.0):
    w = np.exp(lam * np.cos(th))
    W = np.outer(w, w)
    for (a, b, c, d) in [(0, 1, 0, 1), (0, 1, 1, 2), (0, 2, 1, 3)]:
        val = np.trapezoid(np.trapezoid(W * D(a, b) * D(c, d), th, axis=1), th)
        print("  lam=%.1f (%d,%d)x(%d,%d): weighted integral = %+.6e" % (lam, a, b, c, d, val))
