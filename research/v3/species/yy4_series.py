# yy4_series.py -- YY4: exact one-plaquette-graph series of SU(2) plaq, Pade/Borel tests.
# Zero MC. Deterministic. No repo file modified.
from fractions import Fraction as Fr
import math
import numpy as np
from scipy.special import iv
from scipy.integrate import quad

NM = 18

def smul(A, B, n):
    C = [Fr(0)] * (n + 1)
    for i, a in enumerate(A):
        if a:
            for j, b in enumerate(B):
                if b and i + j <= n:
                    C[i + j] += a * b
    return C

def sinv(A, n):
    B = [Fr(0)] * (n + 1)
    B[0] = Fr(1) / A[0]
    for k in range(1, n + 1):
        s = Fr(0)
        for i in range(1, k + 1):
            s += A[i] * B[k - i]
        B[k] = -s / A[0]
    return B

# u/beta = 2*N/(z*D), z=beta^2
# D(z)=sum z^k/(4^k k!(k+1)!), N(z)=sum_{k>=1} z^k/(4^k (k-1)!(k+1)!)
D = [Fr(1, 4**k * math.factorial(k) * math.factorial(k + 1)) for k in range(NM + 1)]
N = [Fr(0)] * (NM + 1)
for k in range(1, NM + 1):
    N[k] = Fr(1, 4**k * math.factorial(k - 1) * math.factorial(k + 1))
Nt = [N[k + 1] for k in range(NM)] + [Fr(0)]
f = [2 * x for x in smul(Nt, sinv(D, NM), NM)]      # f_k = coeff of z^k in plaq/beta
a = f                                            # plaq = sum_k f_k beta^(2k+1)

print("== [1] exact coefficients a_k of plaq = sum a_k beta^(2k+1) (one-plaquette graph) ==")
for k in range(9):
    print("  k=%d  a_k = %s = %.12e" % (k, f[k], float(f[k])))
print("  check a0=1/4 :", f[0] == Fr(1, 4), "  a1=-1/96 :", f[1] == Fr(-1, 96))
print("  ratio |a_k/a_{k+1}| -> radius^2 ; z_pole = -j_{1,1}^2 = %.9f" % (-2.404825557695773**2))
for k in (6, 8, 10, 12):
    print("  k=%2d  |a_k/a_{k+1}| = %.6f" % (k, abs(float(f[k] / f[k + 1]))))

def usp(b):
    return float(iv(2, b) / iv(1, b))

def dusp(b, h=1e-6):
    return (usp(b + h) - usp(b - h)) / (2 * h)

b0 = 2.5
plaq_meas = 0.65264        # repo gh2_pool/gh2_ltA (L=6, 24000 cfg), +-2e-5
print("\n== [2] one-plaquette graph vs anchor at beta=2.5 ==")
print("  I2/I1            = %.9f" % usp(b0))
print("  u_sp'(2.5) FD    = %.9f" % dusp(b0))
# exact-derivative cross-check: d/db(I2/I1) = (I1+I3)/(2I1) - (I2/I1)*(I0+I2)/(2I1)
u0 = usp(b0)
dex = (iv(1, b0) + iv(3, b0)) / (2 * iv(1, b0)) - u0 * (iv(0, b0) + iv(2, b0)) / (2 * iv(1, b0))
print("  u_sp'(2.5) exact = %.9f   (|FD-exact| = %.2e)" % (dex, abs(dex - dusp(b0))))
print("  (I2/I1)/beta     = %.9f   vs 1/4 = %.9f" % (usp(b0) / b0, 0.25))
print("  anchor/beta      = %.9f   (rise above 1/4 = %.4f%%)" % (plaq_meas / b0, 100 * (plaq_meas / b0 / 0.25 - 1)))
print("  anchor           = %.5f  -> gap delta = %.6f  (rel %.4f)" % (plaq_meas, plaq_meas - usp(b0), plaq_meas / usp(b0) - 1))
print("  threshold plaq/beta = %.6f ; margin vs u_sp' = %.3fx" % (plaq_meas / b0, (plaq_meas / b0) / dusp(b0)))

# Pade in z=beta^2
def pade_coeffs(c, m, n):
    A = np.zeros((n, n)); rhs = np.zeros(n)
    for k in range(1, n + 1):
        for j in range(1, n + 1):
            idx = m + k - j
            A[k - 1, j - 1] = c[idx] if idx >= 0 else 0.0
        rhs[k - 1] = -c[m + k]
    q = np.linalg.solve(A, rhs) if n else np.array([])
    qq = np.concatenate(([1.0], q))
    p = [sum(qq[j] * c[i - j] for j in range(0, min(i, n) + 1)) for i in range(m + 1)]
    return p, qq

def pev(p, qq, z):
    return sum(p[i] * z**i for i in range(len(p))) / sum(qq[j] * z**j for j in range(len(qq)))

def pdev(p, qq, z):
    P = sum(p[i] * z**i for i in range(len(p))); Q = sum(qq[j] * z**j for j in range(len(qq)))
    dP = sum(i * p[i] * z**(i - 1) for i in range(1, len(p)))
    dQ = sum(j * qq[j] * z**(j - 1) for j in range(1, len(qq)))
    return (dP * Q - P * dQ) / Q**2

cf = [float(x) for x in f]
z0 = b0**2
print("\n== [3] diagonal Pade [m/m] in z=beta^2, evaluated at z=6.25 ==")
print("   m   Pade(6.25)      err vs I2/I1     u'(2.5)[Pade]")
for m in range(1, 8):
    p, qq = pade_coeffs(cf, m, m)
    v = pev(p, qq, z0); dv = pdev(p, qq, z0)
    up = v + 2 * z0 * dv
    print("  %2d   %.9f   %+.3e    %.6f" % (m, v, v - usp(b0) / b0, up))

# Borel-Pade: B(s)=sum a_k s^(2k+1)/(2k+1)! ; u(beta)=int_0^inf e^-t B(beta t) dt
print("\n== [4] Borel-Pade (g(y)=sum a_k y^k/(2k+1)!, y=s^2) ==")
g = [float(f[k]) / math.factorial(2 * k + 1) for k in range(NM + 1)]
for m in (4, 6, 8):
    p, qq = pade_coeffs(g, m, m)
    integ = lambda t, z=z0: math.exp(-t) * t * pev(p, qq, z * t * t)
    val, err = quad(integ, 0, 60, limit=400)
    print("  [%d/%d]  beta*int = %.9f   (target %.9f, diff %+.2e)" % (m, m, b0 * val, usp(b0), b0 * val - usp(b0)))

print("\n== [5] anchor-matched single-power corrections: what derivative do they force? ==")
delta = plaq_meas - usp(b0)
print("  delta = %+.6f  (anchor - one-plaquette graph)" % delta)
for pp in (5, 7, 9, 13):
    dcont = pp * delta / b0**pp * b0**(pp - 1)   # = pp*delta/beta at beta=2.5
    print("  p=%2d:  u' contribution %+.6f -> u'(2.5) = %+.6f   %s threshold %.6f"
          % (pp, dcont, dusp(b0) + dcont, "ABOVE" if dusp(b0) + dcont > plaq_meas / b0 else "below", plaq_meas / b0))

print("\n== [6] exact arithmetic: can plaq/beta be monotone on [0,2.5]? ==")
print("  plaq/beta at 0+ = 1/4 = %.9f" % 0.25)
print("  plaq/beta at 2.5     = %.9f  (ratio %.6f)" % (plaq_meas / b0, (plaq_meas / b0) / 0.25))
zs = np.linspace(0, 6.25, 200001)
# max of u/beta for the one-plaquette graph
ff = np.array([float(sum(f[k] * (z)**k for k in range(NM + 1))) for z in zs])
imax = int(np.argmax(ff))
print("  one-plaquette u/beta: max %.6f at beta=%.3f ; value at 2.5 %.6f" % (ff[imax], math.sqrt(zs[imax]), usp(b0) / b0))

print("\n== [7] measured tangent/secant ratio (repo SS4-EP:43, rho(2.5)=0.7212) ==")
rho = 0.7212
up_meas = rho * plaq_meas / b0
print("  u'(2.5)_meas = %.6f ; threshold %.6f ; slack %.6f ; margin %.4fx" % (up_meas, plaq_meas / b0, plaq_meas / b0 - up_meas, (plaq_meas / b0) / up_meas))
print("  needed derivative accuracy to keep sign: %.6f (%.1f%% of threshold)" % (plaq_meas / b0 - up_meas, 100 * (plaq_meas / b0 - up_meas) / (plaq_meas / b0)))
