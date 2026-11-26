"""P08 — re-audit beta* = 4.30661 under the SU(2) FUNDAMENTAL-DOMAIN correction found in O24.

O24 (O24-L2-PUSH.md §1 F1) proved that o11's lattice chain drifted out of the SU(2) fundamental
domain: U = exp(i x.sig/2) is a bijection on the ball |x| < 2pi, Haar dU = J(x) dx,
J = 4 sin^2(|x|/2)/|x|^2.  Outside the ball is a repeated cover with no restoring force.

Question here: does the SAME bug infect beta*?
  [1] h1e: is  Var = beta^2 V0  /  W0 = 1/(2Cin2pi)  really PURE ANALYTIC (no link sampling)?
  [2] o07: was R(beta) crossed by EXACT QUADRATURE or by MC -- and on which r-domain?
  [3] three beta* values + the domain-broken control (what the bug would do).
  [4] MC control over beta in [3,5]: fundamental domain vs unbounded chain.

Everything deterministic except [4] (labelled [C]); no fitting of any physical constant.
"""
import math
import numpy as np
from scipy.special import sici

TWO_PI = 2.0 * np.pi
PI = math.pi
EULER = 0.5772156649015328606
L = math.log


def gl(f, a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    t = 0.5 * (b - a) * x + 0.5 * (a + b)
    return 0.5 * (b - a) * float(np.sum(w * f(t)))


# ------------------------------------------------------------------ [1] h1e closed form
print("=" * 96)
print("[1] h1e: is the closed form PURE ANALYTIC?  (slice th=pi, x=rho/2 in (0,pi) <=> rho in (0,2pi))")
Z = gl(lambda x: np.sin(x) ** 2 / x, 0.0, PI, 4000)
N1 = gl(lambda x: np.sin(x) ** 3 / x ** 2, 0.0, PI, 4000)
N2 = gl(lambda x: np.sin(x) ** 4 / x ** 3, 0.0, PI, 4000)
Ci = lambda t: float(sici(t)[1])
Zc = 0.5 * (EULER + L(TWO_PI) - Ci(TWO_PI))                 # = Cin(2pi)/2
N1c = 0.75 * (L(3.0) + Ci(PI) - Ci(3 * PI))
N2c = L(2.0) + Ci(TWO_PI) - Ci(4 * PI)
print(f"    Z  quad = {Z:.15f}   closed 0.5*Cin(2pi)        = {Zc:.15f}   d = {Z-Zc:+.2e}")
print(f"    N1 quad = {N1:.15f}   closed (3/4)[ln3+Ci(pi)-Ci(3pi)] = {N1c:.15f}   d = {N1-N1c:+.2e}")
print(f"    N2 quad = {N2:.15f}   closed ln2+Ci(2pi)-Ci(4pi)   = {N2c:.15f}   d = {N2-N2c:+.2e}")
print(f"    Ci via scipy.special.sici: Ci(pi)={Ci(PI):.12f} Ci(2pi)={Ci(TWO_PI):.12f} "
      f"Ci(3pi)={Ci(3*PI):.12f} Ci(4pi)={Ci(4*PI):.12f}  Cin(2pi)={EULER+L(TWO_PI)-Ci(TWO_PI):.12f}")
V0 = (N2 * Z - N1 ** 2) / (4 * Z ** 2)
V0c = (N2c * Zc - N1c ** 2) / (4 * Zc ** 2)
Cin2 = 2.0 * Zc
W0 = 1.0 / (2.0 * Cin2)
bstar_h1e = math.sqrt(W0 / V0c)
print(f"    V0 (quadrature) = {V0:.15f}    V0 (closed) = {V0c:.15f}    |d| = {abs(V0-V0c):.1e}")
print(f"    W0 = 1/(2Cin2pi) = {W0:.12f}   V0*W0 = {V0c*W0:.9f}  V0/W0 = {V0c/W0:.9f}")
print(f"    beta*(h1e closed) = sqrt(W0/V0) = {bstar_h1e:.8f}")
print("    domain check: the closed form's upper limit x=pi is rho=2pi, i.e. the fundamental-domain")
print("    boundary.  Without that cut the integral DIVERGES (log), so the closed form USES the domain,")
X = [PI, 10.0, 100.0, 1000.0, 1e4]
Ix = [gl(lambda x: np.sin(x) ** 2 / x, 0.0, t, 8000) for t in X]
print("      int_0^X sin^2x/x dx  for X = pi,10,100,1000,1e4 :", " ".join(f"{v:.6f}" for v in Ix))
print(f"      predicted 0.5*Cin(2pi)+0.5*ln(X/pi) at X=1e4      = {Zc+0.5*L(1e4/PI):.6f}")
print("    => finite-by-construction, never sampled: PURE ANALYTIC.  [T]")

# --------------------------------------- [2] o07 crossing, independent exact quadrature
print("=" * 96)
print("[2] o07: R(beta)=Var(F'|0)/E(F''|0) on ONE plaquette, th=pi, slice {x3=0}.")
print("    KERNEL: measure e^{-S} J r dr on (0,2pi); at th=pi S==beta on the slice (Haar).")
print("    o07_ratio.py::n2_exact uses r in (1e-12, 2pi-1e-11), nr=400003  -> EXACT QUADRATURE,")
print("    NOT MC, and the upper limit 2pi IS the fundamental-domain cut.  Re-derived here from")
print("    scratch (midpoint rule, analytic Phi'):")


def slicew(rcap, nr=200001):
    r = (np.arange(nr) + 0.5) * (rcap / nr)          # midpoint rule, no r=0 / r=2pi endpoint
    return r, np.sin(0.5 * r) ** 2 / r                # e^{-S}J r dr  (up to const), S=beta


def moments(beta, rcap=TWO_PI, nr=200001):
    r, w = slicew(rcap, nr)
    w = w / w.sum()
    Phi_r = 2.0 / r - np.cos(0.5 * r) / np.sin(0.5 * r)   # d/dr[-log J]
    F1 = beta * np.sin(0.5 * r) / r                        # <e3,grad S>   (Phi' has no e3 part)
    F2 = Phi_r / r                                         # d^2_v(S+Phi), beta-term = 0 at th=pi
    m1 = (w * F1).sum()
    return (w * F1 * F1).sum() - m1 * m1, (w * F2).sum()


v_lo, e_lo = moments(0.0)
print(f"    beta=0 : Var={v_lo:.3e}  E={e_lo:.12f}   (W0 target {W0:.12f};  d={e_lo-W0:+.1e})")
for b in (3.0, 3.5, 4.0, 4.30661, 4.5, 5.0):
    v, e = moments(b)
    print(f"    beta={b:7.5f}  Var={v:.9f}  E={e:.9f}  R={v/e:.9f}  R-1={v/e-1:+.3e}  Var/b^2={v/b**2:.9f}")


def cross(f, a, b, it=60):
    fa = f(a)
    for _ in range(it):
        m = 0.5 * (a + b)
        if (f(m) > 0) == (fa > 0):
            a, fa = m, f(m)
        else:
            b = m
    return 0.5 * (a + b)


bs_q = cross(lambda b: moments(b)[0] - moments(b)[1], 3.0, 5.0)
v, e = moments(bs_q)
print(f"    beta*(o07 re-derived, fundamental domain r<2pi) = {bs_q:.8f}   "
      f"Var={v:.9f} E={e:.9f} Var-E={v-e:+.2e} R={v/e:.10f}")
# 2D (r,phi) cross-check: is the phi-integral really trivial?
nr, nph = 4001, 64
rr = (np.arange(nr) + 0.5) * (TWO_PI / nr)
pp = (np.arange(nph) + 0.5) * (TWO_PI / nph)
R2, P2 = np.meshgrid(rr, pp, indexing="ij")
Wt = np.sin(0.5 * R2) ** 2 / R2                        # J * r dr dphi  (J=4sin^2/r^2, area=r dr dphi)
Wt = Wt / Wt.sum()
F1 = bs_q * np.sin(0.5 * R2) / R2
F2 = (2.0 / R2 - np.cos(0.5 * R2) / np.sin(0.5 * R2)) / R2
m1 = (Wt * F1).sum()
print(f"    2D (r,phi) check at beta*: Var={((Wt*F1*F1).sum()-m1*m1):.9f}  "
      f"E={(Wt*F2).sum():.9f}   (phi-trivial, 1D reduction exact)")

# ----------------------------------- [3] domain-broken control: beta* vs the escape cutoff
print("=" * 96)
print("[3] CONTROL — what the O24 bug would do to beta*: integrate the SAME kernel to r=R_max")
print("    (repeated-cover measure; total mass ~ 0.5 ln(R_max/2pi), improper as R_max->oo).")
print("      R_max/2pi   int w (mass)   E        Var/b^2     beta*(broken)")
for k in (1, 2, 4, 8, 16, 64, 1024):
    rcap = k * TWO_PI
    r, w = slicew(rcap, 400001)
    mass = w.sum()
    Phi_r = 2.0 / r - np.cos(0.5 * r) / np.sin(0.5 * r)
    b1 = np.sin(0.5 * r) / r
    wn = w / mass
    A = (wn * b1 * b1).sum() - (wn * b1).sum() ** 2
    E = (wn * Phi_r / r).sum()
    print(f"      {k:6d}    {mass:11.3f}   {E:9.6f}  {A:10.7f}   {math.sqrt(E/A):10.4f}")
print("    => WITHOUT the ball cut beta* is cutoff-dependent and diverges (mass ~ ln R_max,")
print("       normalised Var/beta^2 -> 0).  The domain is not a detail: it DEFINES beta*.")

# --------------------------------- [4] MC control beta in [3,5]: domain vs unbounded chain
print("=" * 96)
print("[4] MC control (2D flat slice, random-walk Metropolis, 64 chains).  Exact same kernel.")
print("    V(x1,x2,s) = beta(1 - cos(r'/2)cos(th/2) + sin(r'/2)sin(th/2)u) - log J,  r'=|(x1,x2,s)|")


def V3(x, s, beta, th=PI):
    r = np.sqrt((x * x).sum(1) + s * s)
    u = s / r
    S = beta * (1.0 - np.cos(0.5 * r) * math.cos(0.5 * th)
                + np.sin(0.5 * r) * math.sin(0.5 * th) * u)
    J = 4.0 * np.sin(0.5 * r) ** 2 / r ** 2
    return S - np.log(J)


def logw(x, domain, rcap):
    """2D flat density p(x) = e^{-S} J = 4 sin^2(r/2)/r^2 on the slice (S==beta at th=pi)."""
    r = np.maximum(np.linalg.norm(x, axis=1), 1e-12)
    lw = 2.0 * np.log(np.abs(np.sin(0.5 * r))) - 2.0 * np.log(r)
    if domain:
        return np.where(r < TWO_PI, lw, -np.inf)
    if rcap is not None:
        return np.where(r < rcap, lw, -np.inf)
    return lw


def mc(beta, domain=True, rcap=None, chains=64, sweeps=8000, burn=2000, eps=0.6, seed=8001,
       h=1e-3):
    rng = np.random.default_rng(seed)
    x = rng.normal(0.0, 0.35, (chains, 2))
    cur = logw(x, domain, rcap)
    acc = 0
    F = np.empty((sweeps, chains))
    K = np.empty((sweeps, chains))
    Rs = 0.0
    for t in range(sweeps + burn):
        prop = x + eps * rng.standard_normal(x.shape)
        lp = logw(prop, domain, rcap)
        ok = np.log(rng.random(chains)) < (lp - cur)
        x[ok] = prop[ok]
        cur[ok] = lp[ok]
        acc += ok.mean()
        if t >= burn:
            j = t - burn
            F[j] = (V3(x, h, beta) - V3(x, -h, beta)) / (2 * h)
            K[j] = (V3(x, h, beta) - 2 * V3(x, 0.0, beta) + V3(x, -h, beta)) / h ** 2
            Rs += np.linalg.norm(x, axis=1).mean()
    Ff = F.ravel()                                    # pooled over chains = iid draws
    Kf = K.ravel()
    # jackknife over chains for the errors of Var(F) and E[K]
    n = chains
    jv, je = [], []
    for c in range(n):
        msk = np.ones(n, bool)
        msk[c] = False
        jv.append(F[:, msk].ravel().var())
        je.append(K[:, msk].ravel().mean())
    jv, je = np.array(jv), np.array(je)
    dv = math.sqrt((n - 1) / n * ((jv - jv.mean()) ** 2).sum())
    de = math.sqrt((n - 1) / n * ((je - je.mean()) ** 2).sum())
    return dict(var=Ff.var(), dvar=dv, E=Kf.mean(), dE=de, r=Rs / sweeps,
                rmax=np.linalg.norm(x, axis=1).max(), acc=acc / (sweeps + burn),
                EmF=Ff.mean())


# sanity: which radial weight does the chain actually sample?  (exact expectations)
_r = (np.arange(400001) + 0.5) * (TWO_PI / 400001)
_pr = 2.0 / _r - np.cos(0.5 * _r) / np.sin(0.5 * _r)
for nm, ww in (("sin^2/r (correct)", np.sin(0.5 * _r) ** 2 / _r),
               ("sin^2   (no 1/r)", np.sin(0.5 * _r) ** 2),
               ("sin^2*r", np.sin(0.5 * _r) ** 2 * _r)):
    ww = ww / ww.sum()
    print(f"    exact under w={nm:18s}: <r>={(_r*ww).sum():.4f}  <Phi'/r>={(ww*_pr/_r).sum():.6f}")
print("    (MC samples p_2D = sin^2(r/2)/r^2  <=>  radial w = sin^2(r/2)/r = the first line)")

for b in (3.0, 4.30661, 5.0):
    m = mc(b)
    print(f"    FD   beta={b:7.5f}  Var={m['var']:.6f}+-{m['dvar']:.6f}  E={m['E']:.6f}+-{m['dE']:.6f}"
          f"  R={m['var']/m['E']:.5f}  <r>={m['r']:.3f} rmax={m['rmax']:.3f} acc={m['acc']:.2f}"
          f"  E[F]={m['EmF']:+.4f}  beta*_implied={math.sqrt(m['E']/(m['var']/b**2)):.4f}")
m = mc(4.30661, domain=False)
print(f"    FREE beta= 4.30661  Var={m['var']:.6f}  E={m['E']:.6f}  R={m['var']/m['E']:.5f}  "
      f"<r>={m['r']:.2f} rmax={m['rmax']:.1f} acc={m['acc']:.2f}   (escapes: <r> vs 2pi="
      f"{TWO_PI:.3f})")
m2 = mc(4.30661, domain=False, rcap=8 * TWO_PI)
print(f"    CAP8 beta= 4.30661  Var={m2['var']:.6f}  E={m2['E']:.6f}  R={m2['var']/m2['E']:.5f}  "
      f"<r>={m2['r']:.2f} rmax={m2['rmax']:.2f} acc={m2['acc']:.2f}")

print("=" * 96)
print("SUMMARY")
print(f"  beta*_h1e_closed     = {bstar_h1e:.8f}")
print(f"  beta*_o07_original   = 4.30660785   (o07_ratio.py exact quadrature, r<2pi)")
print(f"  beta*_o07_funddomain = {bs_q:.8f}   (this script, independent quadrature, r<2pi)")
print(f"  |funddomain - original| = {abs(bs_q-4.30660785):.2e}   target 4.30661")
print("  Clay 0%.  No QED.")
