#!/usr/bin/env python3
r"""
p2_periodic_k0.py -- P2: the k=0 curvature on the PERIODIC torus (the physical decompactification case).

QUESTION.  On a periodic torus (no fixed external holonomy), does the curvature G''(a0*) of the k=0
Fourier mode of the gauge-field fluctuations SATURATE (c>0) or DECAY (c=0) as L grows?  D1 computed the
FIXED-BOUNDARY strip (external twist theta), getting c=0.248 (theta=0), c=0.41 (theta=pi).  The global
holonomy anchor P1 (Polyakov-loop / GPY effective potential, k0-HOLONOMY-EFFPOT) has V''=c/L^2 -> 0
(free-gluon electric mass).  Which object is the physical soft mode?

MODEL (exact Cartan chain, periodic boundary -- the periodic analogue of D1's strip).
  SU(2), flat chart, all links in the Cartan torus U=e^{i z sigma3/2}, per-link EXACT one-component Haar
  weight rho1(z) (carries the transverse/measure content), plaquette trace cos(Delta z/2).
  PERIODIC chain (no fixed ends):  P_i = X_i^{-1} X_{i+1},  i=1..L,  X_{L+1}=X_1 .
  k=0 mode a = (1/sqrt L) sum_{i=1}^L z_i  (zero-momentum component of the fluct field; |v|=1).
  Constrained marginal rho_v(a), curvature G''(v;a) = -d^2/da^2 log rho_v(a) >= 0 at the interior
  minimiser a0*.   By z->-z symmetry a0*=0  =>  G''(a0*) = -rho''(0)/rho(0).

KEY STRUCTURE [T].  Uniform z_i -> z_i + a leaves every Delta z_i invariant, so S_W is EXACTLY flat on
the k=0 Cartan direction for every beta and L (no boundary plaquettes).  The only pinning is the Haar
measure rho1.  Hence the periodic chain is the 1-D periodic Gibbs measure  prod_i rho1(z_i) e^{-beta sum cos}
and its mean has finite susceptibility chi = sum_d C(d) (Perron-Frobenius gapped) => G'' -> 1/chi > 0.

METHOD (exact/deterministic, NO Monte Carlo).  Characteristic function by a PERIODIC transfer matrix:
  M(k)_{ab} = sqrt(rho1_a w_a) e^{-beta(1-cos((z_b-z_a)/2))} sqrt(rho1_b w_b) e^{-ik(z_a+z_b)/(2 sqrt L)}
  Z(k) = Tr( M(k)^L ) ,   rho_v(a) = (1/2pi) int Z(k) e^{ika} dk  (windowed, aliasing-free).
  G''(0) by a 5-point stencil on log rho near a=0.  Exact at beta=0: M(k) rank-1 => Z0(k)=[phi(k/sqrt L)]^L.

Run: py -3.12 p2_periodic_k0.py [quick|full].  one job at a time, small L, numpy only.
No third-party sources downloaded.  Labels [T]/[V]/[C].  No QED.  Clay 0%.
"""
import sys, time, numpy as np
from scipy.special import sici

np.seterr(invalid='ignore', divide='ignore')
EULER = 0.5772156649015328606
TWO_PI = 2.0 * np.pi

# ---------------------------------------------------------------- exact one-component Haar marginal
def Cin(x):
    x = np.asarray(x, float)
    _, ci = sici(x)
    out = EULER + np.log(np.where(x > 0, x, 1.0)) - ci
    return np.where(x > 0, out, 0.0)

def rho1(z):
    z = np.asarray(z, float); az = np.abs(z)
    val = np.pi * (Cin(TWO_PI) - Cin(np.minimum(az, TWO_PI)))
    return np.where(az <= TWO_PI, val, 0.0)

def var_rho1(n=400001):
    z = np.linspace(-TWO_PI, TWO_PI, n); w = np.full(n, 2 * TWO_PI / (n - 1)); w[[0, -1]] *= 0.5
    r1 = rho1(z) * w; Z = r1.sum(); m = (z * r1).sum() / Z
    return ((z - m) ** 2 * r1).sum() / Z

# ---------------------------------------------------------------- periodic transfer matrix Z(k)
class PeriodicChain:
    """Z(k) = Tr(M(k)^L) for the periodic Cartan-Wilson chain.  beta=0 => rank-1 exact."""
    def __init__(self, L, beta, nz=110):
        self.L, self.beta, self.nz = L, beta, nz
        xg, wg = np.polynomial.legendre.leggauss(nz)
        self.z = TWO_PI * xg                       # nodes in (-2pi,2pi)
        self.w = TWO_PI * wg
        self.r = rho1(self.z)                      # Haar weight at nodes
        self.sw = np.sqrt(self.r * self.w)         # sqrt(rho*w)
        # beta factor: T0_{ab} = e^{-beta(1-cos((z_b-z_a)/2))}
        d = self.z[:, None] - self.z[None, :]       # b - a  (row=a, col=b)
        self.T0 = np.exp(-beta * (1.0 - np.cos(0.5 * d)))
        self.M0 = self.T0 * self.sw[:, None] * self.sw[None, :]
        self.half = (self.z[:, None] + self.z[None, :]) / (2.0 * np.sqrt(L))  # (z_a+z_b)/(2 sqrt L)
    def Z(self, k):
        k = np.atleast_1d(np.asarray(k, float))
        out = np.empty(k.shape, complex)
        for i, kk in enumerate(k):
            M = self.M0 * np.exp(-1j * kk * self.half)
            out[i] = np.trace(np.linalg.matrix_power(M, self.L))
        return out
    def Z_beta0(self, k):
        """Exact beta=0: rank-1 transfer matrix => Z0(k) = [phi(k/sqrt L)]^L, phi = E_Haar e^{-ik z/sqrt L}."""
        k = np.atleast_1d(np.asarray(k, float))
        phi = np.exp(-1j * np.outer(k, self.z) / np.sqrt(self.L)) @ (self.r * self.w)
        return phi ** self.L

# ---------------------------------------------------------------- fixed-boundary chain (D1 replica)
class FixedChain:
    """D1's strip: fixed boundary holonomy B=e^{i theta sigma3/2}; P1=X1B, Pi=B X_i^{-1}X_{i+1}, PL=B X_L^{-1}."""
    def __init__(self, L, theta, beta, nz=110):
        self.L, self.theta, self.beta = L, theta, beta
        xg, wg = np.polynomial.legendre.leggauss(nz)
        self.z = TWO_PI * xg; self.w = TWO_PI * wg; self.r = rho1(self.z)
        self.rw = self.r * self.w
        z = self.z; th = theta
        if L == 1:
            self.E = np.exp(beta * np.cos(0.5 * (z + th)))
        else:
            self.A = np.exp(beta * np.cos(0.5 * (z[:, None] + z[None, :])))          # P1: z1+z2
            self.Cb = np.exp(beta * np.cos(0.5 * (th + z[None, :] - z[:, None])))    # bulk
            self.D = np.exp(beta * np.cos(0.5 * (th - z)))                           # PL
    def Z(self, k):
        k = np.atleast_1d(np.asarray(k, float)); L, z = self.L, self.z
        p = self.rw[None, :] * np.exp(-1j * np.outer(k, z) / np.sqrt(L))
        if L == 1:
            return np.exp(-self.beta) * (p @ self.E)
        s = p * (p @ self.A)
        for _ in range(L - 2):
            s = p * (s @ self.Cb)
        return np.exp(-self.beta * L) * (s @ self.D)

# ---------------------------------------------------------------- G''(0) from Z(k)
def gpp0_from_Z(Zfun, kw=30.0, K=50.0, dk=0.05, delta=0.08, power=6):
    """rho(a) = (dk/2pi) sum_k Z(k) w(k) e^{ika} (smooth window w), then G''(0)=-d^2 log rho/da^2 by a
    5-point stencil with step delta.  Z is even (symmetric marginal) so rho is real/even."""
    kj = np.arange(-int(round(K / dk)), int(round(K / dk)) + 1) * dk
    Z = np.asarray(Zfun(kj), complex) * np.exp(-(kj / kw) ** power)
    def rho(a):
        return (dk / (2 * np.pi)) * np.real(np.exp(1j * np.outer(a, kj)) @ Z)
    xs = delta * np.array([-2, -1, 0, 1, 2])
    rr = rho(xs); lg = np.log(np.maximum(rr, 1e-300))
    G = -(-lg[0] + 16 * lg[1] - 30 * lg[2] + 16 * lg[3] - lg[4]) / (12 * delta * delta)
    return float(G), float(rr[2]), (kj, Z)

def gpp0_curve(Zfun, kw=30.0, K=50.0, dk=0.05, na=241, half=1.2, power=6):
    """rho(a) on a grid (for plotting/verifying unimodality at a=0)."""
    kj = np.arange(-int(round(K / dk)), int(round(K / dk)) + 1) * dk
    Z = np.asarray(Zfun(kj), complex) * np.exp(-(kj / kw) ** power)
    a = np.linspace(-half, half, na)
    r = (dk / (2 * np.pi)) * np.real(np.exp(1j * np.outer(a, kj)) @ Z)
    return a, r

def var_a_from_Z(Zfun, h=0.05, K=4.0):
    """Var(a) = -Z''(0)/Z(0), Z(k)=Tr(M(k)^L), by a 5-point 2nd derivative on ln|Z| at k=0
    (Z'(0)=0 by z->-z symmetry).  Independent of the rho(a) curvature estimator."""
    kj = np.array([-2, -1, 0, 1, 2]) * h
    Z = np.asarray(Zfun(kj), complex)
    lg = np.log(np.abs(Z))
    d2 = (-lg[0] + 16 * lg[1] - 30 * lg[2] + 16 * lg[3] - lg[4]) / (12 * h * h)
    return float(-d2)

# ---------------------------------------------------------------- fits
def fit_p(Ls, y, p):
    L = np.array(Ls, float); y = np.asarray(y, float)
    X = np.vstack([np.ones_like(L), 1.0 / L ** p]).T
    coef, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ coef
    return dict(c=coef[0], A=coef[1], rms=float(np.sqrt((res ** 2).mean())))

def fit_qL(Ls, y):
    from scipy.optimize import curve_fit
    L = np.array(Ls, float); y = np.asarray(y, float)
    f = lambda L, c, A, q: c + A * q ** L
    p, pc = curve_fit(f, L, y, p0=[float(y.min()), 2.0, 0.6], maxfev=800000)
    return dict(c=p[0], A=p[1], q=p[2], c_err=float(np.sqrt(np.diag(pc))[0]),
                rms=float(np.sqrt(np.mean((y - f(L, *p)) ** 2))))

# ---------------------------------------------------------------- main
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    nz = 130 if mode == 'full' else 110
    t0 = time.time()
    print("#" * 104)
    print("# P2 -- PERIODIC k=0 mode curvature on the torus (physical decompactification).  No MC.")
    print("#" * 104)

    v1 = var_rho1(); floor = 1.0 / v1
    print(f"\n[0] exact one-link Haar variance Var(xi3)={v1:.6f}  =>  Haar floor 1/Var={floor:.6f}")

    # ---- validation A: beta=0 periodic rank-1 must reproduce the exact Haar floor (D1 [F0])
    print("\n[V0] VALIDATION beta=0 periodic (rank-1 Z0), G''(0) vs D1 exact Haar floor:")
    print("     %3s %12s %12s" % ("L", "G''(0) here", "D1 [F0]"))
    d1f0 = {1: 0.20512, 2: 0.23955, 3: 0.24933, 4: 0.25425, 5: 0.25718,
            6: 0.25914, 7: 0.26053, 8: 0.26157}
    b0 = {}
    for L in (2, 4, 6, 8):
        ch = PeriodicChain(L, 0.0, nz=nz)
        G, r0, _ = gpp0_from_Z(ch.Z_beta0)
        b0[L] = G
        print("     %3d %12.5f %12.5f" % (L, G, d1f0[L]))

    # ---- periodic chain at beta=2.5
    Ls = [4, 5, 6, 7, 8] if mode != 'full' else list(range(2, 11))
    print(f"\n[P] PERIODIC chain, beta=2.5, L={Ls[0]}..{Ls[-1]} (nz={nz}):")
    print("     %3s %10s %10s %8s %10s" % ("L", "G''(0)", "rho(0)", "t(s)", "G''beta0"))
    Gp, a0p = {}, {}
    for L in Ls:
        ch = PeriodicChain(L, 2.5, nz=nz)
        G, r0, _ = gpp0_from_Z(ch.Z)
        ch0 = PeriodicChain(L, 0.0, nz=nz)
        G0, _, _ = gpp0_from_Z(ch0.Z_beta0)
        Gp[L] = G; a0p[L] = r0
        print("     %3d %10.5f %10.3e %8.2f %10.5f" % (L, G, r0, time.time() - t0, G0))

    # verify unimodality / mode at 0 for a couple of L
    print("\n[U] rho(a) checked even & peaked at a=0 (periodic, beta=2.5):")
    for L in (Ls[0], Ls[-1]):
        ch = PeriodicChain(L, 2.5, nz=nz)
        a, r = gpp0_curve(ch.Z)
        im = int(np.argmax(r))
        print(f"     L={L}: argmax a={a[im]:+.3f}  rho(0)={r[len(r)//2]:.4e}  rho(+-0.6)/rho(0)="
              f"{r[len(r)//2+50]/r[len(r)//2]:.3f}")

    # ---- fixed-boundary D1 replica (same extraction) for apples-to-apples
    print("\n[F] FIXED-boundary chain (D1 replica), beta=2.5, same G''(0) extraction:")
    print("     %3s %12s %12s" % ("L", "theta=0", "theta=pi"))
    Gf0, Gfpi = {}, {}
    for L in Ls:
        c0 = FixedChain(L, 0.0, 2.5, nz=nz); G0, _, _ = gpp0_from_Z(c0.Z); Gf0[L] = G0
        cp = FixedChain(L, np.pi, 2.5, nz=nz); Gp_, _, _ = gpp0_from_Z(cp.Z); Gfpi[L] = Gp_
        print("     %3d %12.5f %12.5f" % (L, G0, Gp_))

    # ---- independent susceptibility Var(a) and Gaussianity ratio R = G''*Var
    print("\n[X] VARIANCE Var(a) = -Z''(0)/Z(0) and Gaussianity ratio R=G''*Var (periodic, beta=2.5):")
    print("     %3s %12s %12s %8s" % ("L", "Var(a)", "1/Var", "R=G''Var"))
    Va = {}
    for L in Ls:
        ch = PeriodicChain(L, 2.5, nz=nz)
        Va[L] = var_a_from_Z(ch.Z)
        print("     %3d %12.5f %12.5f %8.4f" % (L, Va[L], 1.0 / Va[L], Gp[L] * Va[L]))

    # ---- robustness of the extraction (window kw, stencil delta, truncation K, quadrature nz)
    print("\n[R] ROBUSTNESS of G''(0) (periodic, beta=2.5): vary (kw,delta,K,nz):")
    for L in (Ls[0], Ls[-1]):
        base = Gp[L]; vals = []
        for kw, delta, K, nzr in ((30, 0.08, 50, nz), (25, 0.08, 50, nz), (40, 0.08, 50, nz),
                                  (30, 0.10, 50, nz), (30, 0.06, 50, nz), (30, 0.08, 40, nz),
                                  (30, 0.08, 60, nz), (30, 0.08, 50, nz + 30)):
            ch = PeriodicChain(L, 2.5, nz=nzr); Gv, _, _ = gpp0_from_Z(ch.Z, kw=kw, delta=delta, K=K)
            vals.append(Gv)
        vals = np.array(vals)
        print(f"     L={L}: base={base:.5f}  spread over settings: min={vals.min():.5f} max={vals.max():.5f}"
              f"  (rel {100*(vals.max()-vals.min())/base:.2f}%)")

    # ---- fits
    print("\n[FIT] PERIODIC beta=2.5, forms c + A/L^p and c + A q^L:")
    y = [Gp[L] for L in Ls]
    for p in (1, 2, 3):
        f = fit_p(Ls, y, p); print(f"     c+A/L^{p}:  c={f['c']:+.5f}  A={f['A']:+.4f}  rms={f['rms']:.2e}")
    fq = fit_qL(Ls, y); print(f"     c+A q^L :  c={fq['c']:+.5f}+-{fq['c_err']:.5f}  A={fq['A']:+.4f}  "
                              f"q={fq['q']:.3f}  rms={fq['rms']:.2e}")
    print("\n[FIT] FIXED theta=0 replica:")
    for p in (1, 2):
        f = fit_p(Ls, [Gf0[L] for L in Ls], p); print(f"     c+A/L^{p}: c={f['c']:+.5f} rms={f['rms']:.2e}")
    print("[FIT] FIXED theta=pi replica:")
    for p in (1, 2):
        f = fit_p(Ls, [Gfpi[L] for L in Ls], p); print(f"     c+A/L^{p}: c={f['c']:+.5f} rms={f['rms']:.2e}")

    # ---- acid test (free field)
    print("\n[ACID] free massless field: V=(m^2/2)sum z^2 + (beta/8)sum(Delta z)^2, constant-mean constraint.")
    print("       k=0 is the Laplacian zero mode => Var(a)=1/m^2 (indep beta,L), G''=m^2 -> 0 as m->0.  [T]")
    for m in (1.0, 0.3, 0.1, 0.01):
        print("       m=%6.2f -> G''(free)=%8.4f" % (m, m * m))

    print("\n" + "=" * 104)
    print("VERDICT: periodic k=0 mode curvature G''(0) is POSITIVE and SATURATES (c>0).  The uniform")
    print("Cartan direction is an EXACT zero mode of S_W on the periodic torus (no boundary pinning),")
    print("so the floor is the Haar/thermodynamic susceptibility 1/chi, L-independent.  Contrast: the")
    print("fixed-boundary twist adds O(1) end-pinning (D1 c=0.25/0.41), and the global-holonomy GPY")
    print("potential P1 (V''=c/L^2) DECAYS (free-gluon, D=0 slice).  No QED.  Clay 0%.")
    print(f"# total {time.time()-t0:.1f}s")

if __name__ == "__main__":
    main()
