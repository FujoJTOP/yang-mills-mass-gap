import numpy as np
from scipy.special import iv

# FIXED normalization: use omega in [0,pi] on S^3 (NOT theta/2).
#   q = (1/2)tr U = cos(omega)
#   chi_j(omega) = sin((2j+1)omega)/sin(omega)      (= U_{2j}(cos omega))
#   Haar: int f dU = (2/pi) int_0^pi f(omega) sin^2(omega) domega
#   f = sum_j c_j chi_j  with  c_j = int f chi_j dU

def charval(j2, om):
    s = np.sin(om)
    num = np.sin((j2 + 1.0) * om)
    return np.where(np.abs(s) < 1e-12, (j2 + 1.0), num / np.where(np.abs(s) < 1e-12, 1.0, s))

def ccoeffs(betap, jmax2, nq=12001):
    om = np.linspace(0, np.pi, nq)
    q = np.cos(om)
    f = np.exp(-betap * (1.0 - q))
    meas = (2.0 / np.pi) * np.sin(om) ** 2
    out = []
    for j2 in range(0, jmax2 + 1):
        out.append(np.trapezoid(f * charval(j2, om) * meas, om))
    return np.array(out)

print("=== c_0 must equal e^{-b'} 2I_1(b')/b' ===")
for betap in (0.5, 1.5, 3.0):
    c = ccoeffs(betap, 4)
    exact = np.exp(-betap) * 2 * iv(1, betap) / betap
    print("  betap=%.1f : c_0(quad)=%.10f   e^{-b}2I1(b)/b=%.10f   rel.err=%.2e"
          % (betap, c[0], exact, abs(c[0] - exact) / exact))

print()
print("=== reconstruction check ===")
betap = 1.5
c = ccoeffs(betap, 30)
om = np.array([0.0, 0.7, 1.4, 2.1, 2.8])
ft = np.exp(-betap * (1.0 - np.cos(om)))
fs = sum(c[j2] * charval(j2, om) for j2 in range(len(c)))
print("  omega :", om)
print("  f_true:", np.array2string(ft, precision=8))
print("  f_ser :", np.array2string(fs, precision=8))
print("  max rel err = %.2e" % np.max(np.abs(fs - ft) / ft))

print()
print("=== positivity of c_j (f > 0 => all c_j > 0) ===")
for betap in (0.5, 1.5, 3.0):
    c = ccoeffs(betap, 12)
    print("  betap=%.1f : c_j(j2=0..6) = %s   min = %+.3e"
          % (betap, np.array2string(c[:7], precision=6), c[:7].min()))
