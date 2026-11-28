import numpy as np
from scipy.special import iv

# Settle the closed form for the character coefficients of f = exp(-b(1-cos w)).
#   c_j = (2/pi) int_0^pi f(cos w) * sin((2j+1)w)/sin(w) * sin^2(w) dw
#       = (2/pi) int_0^pi f(cos w) sin((2j+1)w) sin(w) dw
# Claim:  c_j = e^{-b} ( I_{2j}(b) - I_{2j+2}(b) )
#              = e^{-b} * ( 2(2j+1)/b ) * I_{2j+1}(b)        [Bessel recurrence]

def charval(j2, om):
    s = np.sin(om); num = np.sin((j2 + 1.0) * om)
    return np.where(np.abs(s) < 1e-12, (j2 + 1.0), num / np.where(np.abs(s) < 1e-12, 1.0, s))

def cquad(b, jmax2, nq=200001):
    om = np.linspace(0, np.pi, nq)
    f = np.exp(-b * (1.0 - np.cos(om)))
    meas = (2.0 / np.pi) * np.sin(om) ** 2
    return np.array([np.trapezoid(f * charval(j2, om) * meas, om) for j2 in range(jmax2 + 1)])

print("=== Bessel values at b=1.5 (recompute from scratch) ===")
b = 1.5
for n in range(6):
    print("  I_%d(%.1f) = %.10f" % (n, b, iv(n, b)))
print("  I_0 - I_2      = %.10f" % (iv(0, b) - iv(2, b)))
print("  (2/1.5)*I_1    = %.10f" % (2 / b * iv(1, b)))
print("  => recurrence I_0 - I_2 = (2/b) I_1 : %s"
      % ("OK" if abs((iv(0,b)-iv(2,b)) - 2/b*iv(1,b)) < 1e-12 else "MISMATCH"))
print()

print("=== compare quadrature c_j with the closed form ===")
print("%5s %6s %16s %16s %16s %10s" % ("j2", "b", "c_j (quad)", "e^-b(I_j-I_j+2)", "e^-b 2(j+1)/b I_j+1", "max rel"))
for bb in (0.5, 1.5, 3.0):
    c = cquad(bb, 8)
    for j2 in range(0, 9):
        j = j2 / 2.0
        f1 = np.exp(-bb) * (iv(2 * j, bb) - iv(2 * j + 2, bb))
        f2 = np.exp(-bb) * (2 * (2 * j + 1) / bb) * iv(2 * j + 1, bb)
        rel = max(abs(c[j2] - f1) / abs(f1), abs(f1 - f2) / abs(f2))
        print("%5d %6.2f %16.10f %16.10f %16.10f %10.2e"
              % (j2, bb, c[j2], f1, f2, rel))
    print()
