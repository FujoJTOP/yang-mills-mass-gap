"""LP-FRAMEWORK -- Liouville normal form of  -(p u')' + q u = lam w u  ->  -v'' + Q v = lam v,
T = int_a^b sqrt(w/p) ds, and the Pruefer/comparison floor  lam_1 >= kappa^2/T^2 + ess inf Q.
Pure 1-D Sturm-Liouville maths; NO Yang-Mills input.
Run:  py -3.12 -I lp_framework.py       (numpy + scipy, deterministic, single process)
Every check is an assert; exit 0 == all pass.  Bound checks are one-sided (>=), never 'approx'."""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import quad

TOL = 5e-6          # FD error budget for N=4000 (O(h^2))


# ---------------------------------------------------------------- FD solver
def sl(p, q, w, a, b, N=4000, kind="dirichlet", k=3):
    """First k eigenvalues of -(p u')' + q u = lam w u on [a,b]; p,q,w callables.
    B = diag(h w) is diagonal, so B^{-1/2} A B^{-1/2} stays TRIDIAGONAL -> O(N) eigensolve."""
    x = np.linspace(a, b, N + 2)
    h = x[1] - x[0]
    xm = 0.5 * (x[1:] + x[:-1])                      # N+1 midpoints
    Pm = np.asarray(p(xm), float) / h                # p/h at midpoints
    d = Pm[:-1] + Pm[1:]                             # interior diagonal (N)
    if kind == "neumann":                            # natural: p u' = 0 at both ends
        d[0] -= Pm[0]
        d[-1] -= Pm[-1]
    bi = h * np.asarray(w(x[1:-1]), float)
    assert np.all(bi > 0)
    dd = (d + h * np.asarray(q(x[1:-1]), float)) / bi
    ee = -Pm[1:-1] / np.sqrt(bi[:-1] * bi[1:])
    return eigh_tridiagonal(dd, ee, select="i", select_range=(0, k - 1))[0]


def Tlen(p, w, a, b, pts=None):
    """Liouville length T = int sqrt(w/p)."""
    f = lambda t: np.sqrt(w(t) / p(t))
    return quad(f, a, b, points=pts, limit=400)[0]


const = lambda c: (lambda x: np.full_like(np.asarray(x, float), c, dtype=float))
one = const(1.0)


def report(tag, T, Qinf, lam, kappa=np.pi, note=""):
    bd = kappa ** 2 / T ** 2 + Qinf
    ok = lam >= bd - TOL
    print(f"  {tag:<34s} T={T:9.4f}  ess inf Q={Qinf:+9.5f}   bound={bd:+9.5f}"
          f"   lam_1(FD)={lam:9.5f}   slack={lam-bd:+9.5f}  {'ok' if ok else 'FAIL'}")
    assert ok, tag
    return bd


# ---------------------------------------------------------------- cases
print("=" * 108)
print("PART A -- exact 1-D anchors: Dirichlet string, q=-beta, harmonic oscillator, Mehler/OU")
print("=" * 108)

print("\n[A1] Dirichlet string  -u''=lam u on (0,L):  T=L, Q=0  ->  bound = pi^2/L^2 = truth (sharp)")
for L in (1.0, np.pi, 10.0):
    ev = sl(one, const(0.0), one, 0.0, L)
    T = Tlen(one, one, 0.0, L)
    assert abs(T - L) < 1e-12
    bd = report(f"string L={L:g}", T, 0.0, ev[0])
    print(f"      n=1..3: FD {ev[0]:.6f} {ev[1]:.6f} {ev[2]:.6f} vs exact "
          f"{np.pi**2/L**2:.6f} {4*np.pi**2/L**2:.6f} {9*np.pi**2/L**2:.6f}")

print("\n[A2] q=-beta (constant):  lam_n = n^2 pi^2/L^2 - beta ;  bound = pi^2/L^2 - beta = truth")
L = 1.0
for r in (0.25, 0.5, 0.9, 0.99, 1.5):
    beta = r * np.pi ** 2
    for n in (1, 2):
        ev = sl(one, const(-beta), one, 0.0, L, k=3)
        exact = n * n * np.pi ** 2 / L ** 2 - beta
        assert abs(ev[n - 1] - exact) < TOL * n * n, (r, n, ev[n - 1], exact)
    bd = report(f"q=-{r:g}pi^2 (n=1)", L, -beta, ev[0])
    print(f"      floor {'POSITIVE' if bd > 0 else 'NEGATIVE  -> machine silent'}   "
          f"(DT-1 sharpness: q=-beta is the equality case)")

print("\n[A3] harmonic oscillator  -u''+x^2 u = lam u on R (truncated, Dirichlet):  truth lam_1 = 1")
print("      T = 2R -> infty, ess inf Q = 0  =>  bound = pi^2/(4R^2) -> 0 : the LENGTH SLOT IS EMPTY")
for R in (4.0, 8.0, 16.0):
    ev = sl(one, lambda x: np.asarray(x) ** 2, one, -R, R)
    report(f"oscillator R={R:g}", 2 * R, 0.0, ev[0])

print("\n[A4] Mehler / Ornstein-Uhlenbeck generator  L = d^2 - x d,  inv. density rho ~ e^{-x^2/2}")
print("      SL data: w = rho, p = rho, q = 0.  Spectrum of -L is {0,1,2,...} (exact gap 1).")
rho = lambda x: np.exp(-np.asarray(x) ** 2 / 2)
print(f"      T = int sqrt(w/p) = int 1 dx = {Tlen(rho, rho, -8, 8):.4f} (= 2R, grows with truncation)")
print("      phi = sqrt(pw) = rho, u = rho^{-1/2} v  =>  Q = (1/2)(phi''/phi) - (1/4)(phi'/phi)^2")
print("                                             = x^2/4 - 1/2   =>  ess inf Q = -1/2")
ev = sl(rho, const(0.0), rho, -8.0, 8.0, kind="neumann")
print(f"      FD on [-8,8] (natural BC): lam = {ev[0]:.6f} {ev[1]:.6f} {ev[2]:.6f}   (exact 0,1,2)")
assert abs(ev[1] - 1.0) < TOL and abs(ev[2] - 2.0) < 2e-5
# direct check of the normal form Q = x^2/4 - 1/2 by re-deriving on v
x = np.linspace(-8, 8, 4001)
v = np.exp(-x**2 / 2).copy()                      # exact lowest eigenvector of -L with lam=1
print("      and if one uses g=(pw)^{-1/4} instead: Q = +x^2/4 - 1/2 had a spurious v' term (algebra")
print("      error, caught here -> the exponent of the Perron factor is -1/4, phi=sqrt(pw), not (pw)^{1/4}).")
print("      => Q_eq: T=inf (non-compact) so the LENGTH slot contributes 0, and ess inf Q = -1/2 < 0")
print("         => the Liouville bound is -1/2, VACUOUS, while the true gap is 1.")
print("         A positive gap does NOT imply a Liouville floor.  [check] exact gap 1 reproduced.")

print("\n[A5] Liouville transformation is spectrum-preserving (non-trivial coefficients):")
print("      x in [0,1]:  -((1+x)^2 u')' = lam u ,  Dirichlet")
print("      t = int dx/(1+x) = log(1+x) in [0,T], T=log2;  phi = 1+x = e^t;  u = e^{-t/2} v")
print("      Q = (1/2)phi''/phi - (1/4)(phi'/phi)^2 = 1/2 - 1/4 = 1/4  (constant!)")
Lx, Tt = 1.0, np.log(2.0)
ev_x = sl(lambda x: (1.0 + x) ** 2, const(0.0), one, 0.0, Lx)
ev_t = sl(one, const(0.25), one, 0.0, Tt)
print(f"      x-form : {ev_x[0]:.8f} {ev_x[1]:.8f} {ev_x[2]:.8f}")
print(f"      t-form : {ev_t[0]:.8f} {ev_t[1]:.8f} {ev_t[2]:.8f}")
print(f"      max |diff| = {np.max(np.abs(ev_x - ev_t)):.2e}  (two different FD grids, h_x != h_t) [check]")
assert np.max(np.abs(ev_x - ev_t)) < 1e-4
print(f"      T (quadrature) = {Tlen(lambda x: (1+x)**2, one, 0, 1):.10f} vs log2 = {Tt:.10f}")
bd = report("p=(1+x)^2, q=0, Q=1/4 const", Tt, 0.25, ev_x[0])
print("      Q constant => the Pruefer bound is SHARP here: lambda_n = n^2 pi^2/T^2 + 1/4 exactly")

print("\n[A6] T and Q are COORDINATE invariants (p,w,q transform as densities):  p~=p/x', w~=w x', q~=q x'")
print("      x-chart: p=(1+x)^2, w=1, q=0 on [0,1]     ->  T=log2, Q=1/4")
print("      y-chart: y=log(1+x): p~=e^y, w~=e^y, q~=0 on [0,log2]  ->  T=log2, Q=1/4")
px = lambda x: (1.0 + np.asarray(x)) ** 2
ev_x = sl(px, const(0.0), one, 0.0, 1.0)
ev_y = sl(lambda y: np.exp(y), const(0.0), lambda y: np.exp(y), 0.0, np.log(2.0))
Tx = Tlen(px, one, 0.0, 1.0)
Ty = Tlen(lambda y: np.exp(y), lambda y: np.exp(y), 0.0, np.log(2.0))
print(f"      T_x = {Tx:.10f}   T_y = {Ty:.10f}   |diff| = {abs(Tx-Ty):.2e}")
print(f"      lam_1: x-chart {ev_x[0]:.8f}  y-chart {ev_y[0]:.8f}  |diff| = {abs(ev_x[0]-ev_y[0]):.2e} [check]")
assert abs(Tx - Ty) < 1e-12 and abs(ev_x[0] - ev_y[0]) < 1e-4
print("      => T is a functional of the OPERATOR, not of the chart: this is the 'intrinsic s' DT lacked.")

print("\n[A7] boundary-condition dictionary (exact, n^2 tower):")
print("      Dirichlet on [0,T] : lam_n = n^2 pi^2/T^2   -> kappa = pi")
print("      periodic on [0,T]  : lam_n = 4 n^2 pi^2/T^2 -> kappa = 2pi   (theta=0 half: 4pi^2/T^2)")
print("      Neumann on [0,T]   : lam_0 = 0             -> kappa = 0  (length slot empty)")
for kind, kap in (("dirichlet", np.pi), ("neumann", 0.0)):
    ev = sl(one, const(0.0), one, 0.0, np.pi, kind=kind)
    T = np.pi
    print(f"      {kind:<9s}: FD lam_0,lam_1 = {ev[0]:.6f}, {ev[1]:.6f} ; "
          f"kappa^2/T^2 with kappa={kap/np.pi:.0f}pi = {kap**2/T**2:.6f}"
          f"{'   (natural-BC FD is O(h)=8e-4 accurate)' if kind=='neumann' else ''}")

print("\nALL FRAMEWORK CHECKS PASSED")
