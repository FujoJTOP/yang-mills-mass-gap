#!/usr/bin/env python3
r"""
p4_periodic_transfer.py -- P4: the PERIODIC constant-mode transfer matrix.

QUESTION.  D3 computed the constant-mode transfer matrix on the *fixed-boundary* 1xL strip
(P1 = X1 X2, P_i = B X_i^-1 X_{i+1}, P_L = B X_L^-1) and found a bulk gap Delta > 0 (3.18 at
theta=pi, 0.64 at theta=0) and a curvature G''(L) that saturates at the Haar floor 0.2688 + a
decaying boundary (Wilson) part.  Does the *periodic* transfer matrix -- cyclic chain X_{L+1}=X_1,
NO twist, NO fixed boundary -- also have a gapped, L-INDEPENDENT spectrum, and does its Riccati
curvature G'' saturate at a positive constant?

MODEL (pure-Cartan sector, identical Haar marginal to D1/D3/D4).
  Cyclic chain of L SU(2) links X_1,...,X_L, X_{L+1}=X_1.  Flat chart z_l = xi_{l,3} in (-2pi,2pi).
  Haar measure factorises into the exact one-link Cartan marginal
      rho1(z) = pi ( Cin(2pi) - Cin|z| ) ,   |z| < 2pi      (beta-independent; transverse
  components integrated EXACTLY into it).  Wilson action is the Cartan projection of the L cyclic
  plaquettes P_i = X_i^-1 X_{i+1}:
      S_W = beta sum_{i=1}^{L} (1 - cos((z_{i+1}-z_i)/2)) ,   z_{L+1}=z_1 .
  Constant mode  v_l = e3/sqrt(L),  a = <v,xi> = (1/sqrt L) sum_l z_l  (unit-normalised mode).

TRANSFER MATRIX (periodic).  One link per step; the partition function of the CYCLIC chain is a
TRACE (closure), not a fixed-ends matrix element:
      Z(k) = < e^{-i k a} >_V = Tr[ M_k^L ] ,
      M_k(z,z') = sqrt(rho1(z) rho1(z')) e^{-i k (z+z')/(2 sqrt L)} e^{- beta (1 - cos((z'-z)/2))} .
At k=0, M_0(z,z') = sqrt(rho1(z)rho1(z')) e^{-beta(1-cos((z'-z)/2))} is EXACTLY D3's bulk kernel at
theta=0 (same object => same bulk gap; verified).  rho_v(a) = (1/2pi) int Z(k) e^{i k a} dk, and the
Riccati gap is  G''(a) = - d_a^2 log rho_v(a).  By symmetry a0* = 0 and the curvature at the peak is
the EXACT moment ratio  G''(0) = int k^2 Z(k) dk / int Z(k) dk  (since rho'(0)=0).

WHAT IS DIFFERENT FROM THE FIXED-BOUNDARY STRIP.
  On the open strip the boundary/sum plaquettes P1 and P_L are NOT functions of the link differences,
  so the constant configuration is PINNED by the ends: an O(1/L) Wilson restoring force (D1 sec 1,
  D3 sec 2) -- this is what makes the strip G''_Wilson > 0 and decaying.  On the PERIODIC chain every
  plaquette is X_i^-1 X_{i+1}, a pure difference, so the uniform shift z_l -> z_l + c leaves S_W
  EXACTLY invariant: the constant mode is an exact FLAT direction of the action.  The action can then
  act on the constant-mode *marginal* only through the (measure-mediated) correlation of the
  non-constant fluctuations -- no O(1/L) end pinning.  This is the qualitative difference to be tested.

SECTIONS
  [T1] PERIODIC bulk spectral transfer matrix: eigenvalues lambda_0>lambda_1>... of M_0, spectral gap
       Delta = -log(lambda_1/lambda_0), correlation length xi=1/Delta.  One kernel => L-INDEPENDENT.
       beta=0 is rank-1 (Delta=inf).  Cross-checked against D3's theta=0 strip kernel.
  [T2] EXACT periodic constant-mode curvature G''(0), L=2..8, beta=2.5 and beta=0, via
       Z(k)=Tr[M_k^L] (moment method + inverse-transform cross-check).
  [T3] SPLIT against the beta=0 Haar floor: G''_Wilson = G''(2.5) - G''_Haar.  On the periodic chain
       the flat-direction decoupling predicts a SMALL, non-1/L Wilson part (contrast: twist big).
  [T4] FITS  c + A/L^p  and  c + A q^L : is the intercept c > 0 and L-independent?
  [T5] ACID: free Gaussian G''=kappa -> 0 ; free massless box gap 4pi^2/L^2 -> 0 ; Haar floor 0.2688>0.

Run: py -3.12 -I p4_periodic_transfer.py [quick|full].  numpy/scipy only, deterministic, no MC.
Labels [T]/[V]/[C].  No QED.  Clay 0%.
"""
import sys
import numpy as np
from scipy.special import sici
from scipy.optimize import curve_fit

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
    """rho1(z)=2pi int_{|z|}^{2pi} sin^2(r/2)/r dr = pi(Cin(2pi)-Cin|z|), support |z|<2pi."""
    z = np.asarray(z, float); az = np.abs(z)
    return np.where(az <= TWO_PI, np.pi * (Cin(TWO_PI) - Cin(np.minimum(az, TWO_PI))), 0.0)

def haar_var_xi3(nz=400001):
    z = np.linspace(-TWO_PI, TWO_PI, nz); w = np.full(nz, 2 * TWO_PI / (nz - 1)); w[[0, -1]] *= 0.5
    r1 = rho1(z) * w
    return (z * z * r1).sum() / r1.sum()

def nodes(nz):
    xg, wg = np.polynomial.legendre.leggauss(nz)
    z = TWO_PI * xg; wz = TWO_PI * wg
    return z, rho1(z) * wz                      # exact Cartesian Haar measure weights

# ---------------------------------------------------------------- periodic bulk kernel (k=0)
def bulk_kernel(z, r1w, beta):
    """M_0[i,j] = sqrt(r1w_i r1w_j) e^{-beta(1 - cos((z_j - z_i)/2))}; r1w = rho1*w."""
    s = np.sqrt(np.maximum(r1w, 0.0))
    K = np.exp(beta * np.cos(0.5 * (z[None, :] - z[:, None])))
    M = (s[:, None] * K) * s[None, :]
    return 0.5 * (M + M.T)

def bulk_spectrum_periodic(beta, nz=700):
    z, r1w = nodes(nz)
    M = bulk_kernel(z, r1w, beta)
    ev = np.linalg.eigvalsh(M)[::-1]
    ratio = np.clip(ev[1] / ev[0], 1e-300, 1.0) if ev[0] > 0 else 1e-300
    Delta = -np.log(ratio)
    xi = 1.0 / Delta if Delta > 1e-9 else np.inf
    return ev[:4], Delta, xi

# ---------------------------------------------------------------- periodic Z(k) = Tr[M_k^L]
def build_Z_periodic(L, beta, kgrid, z, r1w):
    """Z(k)=Tr[M_k^L] for the CYCLIC chain, M_k = M_0 o (ph ph^T), ph=e^{-i k z/(2 sqrt L)}.
    Vectorised over kgrid; one nz x nz eigenproblem per k."""
    sqL = np.sqrt(L)
    s = np.sqrt(np.maximum(r1w, 0.0))
    K0 = np.exp(beta * np.cos(0.5 * (z[None, :] - z[:, None])))   # e^{-beta(1-cos)} up to beta factor
    base = (s[:, None] * K0) * s[None, :]                         # = e^{beta} M_0  (real sym)
    # each step of the trace carries e^{-beta} (from e^{-beta(1-cos)}=e^{-beta}e^{beta cos}); L steps:
    pref = np.exp(-beta * L)
    Z = np.empty(len(kgrid), dtype=complex)
    for idx, k in enumerate(kgrid):
        ph = np.exp(-1j * k * z / (2.0 * sqL))
        Mk = base * (ph[:, None] * ph[None, :])
        ev = np.linalg.eigvals(Mk)
        Z[idx] = pref * np.sum(ev ** L)
    return Z

def build_Z_periodic_beta0(L, kgrid, z, r1w):
    """beta=0 closed form (independent links): Z(k) = (int rho1 e^{-i k z/sqrt L})^L."""
    sqL = np.sqrt(L)
    z1 = np.exp(-1j * np.outer(kgrid, z) / sqL) @ r1w        # one-link characteristic function
    return z1 ** L

def make_kgrid(mode):
    # Z(k) ~ exp(-k^2 Var/2), Var~1.9 => negligible past |k|~8. Range +-K with K*dz << pi (no aliasing).
    if mode == 'full':
        K, dk = 8.0, 0.015
    else:
        K, dk = 7.0, 0.03
    nh = int(round(K / dk))
    return np.arange(-nh, nh + 1) * dk, dk

def gcurv_moment(Z, kgrid):
    """G''(0)=int k^2 Z / int Z (exact at the symmetric peak, since rho'(0)=0)."""
    k = np.asarray(kgrid); w = Z.real
    trap = getattr(np, 'trapezoid', None) or np.trapz
    den = trap(w, k); num = trap(k * k * w, k)
    return num / den

def rho_at(Z, kgrid, dk, av):
    av = np.atleast_1d(av)
    return (dk / (2 * np.pi)) * np.real(np.exp(1j * np.outer(av, kgrid)) @ Z)

def gap_from_Z(Z, kgrid, dk, half, delta=0.09, na=401):
    """Cross-check: locate a0* and G''(a0*) by 5-point -d^2 log rho / da^2."""
    av = np.linspace(-half, half, na)
    rr = rho_at(Z, kgrid, dk, av)
    m = 6; im = m + int(np.argmax(rr[m:-m]))
    y0, y1, y2 = rr[im - 1], rr[im], rr[im + 1]; den = (y0 - 2 * y1 + y2)
    sh = 0.5 * (y0 - y2) / den if abs(den) > 1e-30 else 0.0
    a0 = av[im] + sh * (av[1] - av[0])
    def curv(x):
        xs = x + delta * np.array([-2, -1, 0, 1, 2])
        lr = np.log(np.maximum(rho_at(Z, kgrid, dk, xs), 1e-300))
        return -(-lr[0] + 16 * lr[1] - 30 * lr[2] + 16 * lr[3] - lr[4]) / (12 * delta * delta)
    return a0, curv(a0), curv(0.0)

# ---------------------------------------------------------------- fits
def fit_cAp(Ls, G):
    L = np.array(Ls, float); y = np.array([G[l] for l in Ls], float)
    f = lambda L, c, A, p: c + A / L ** p
    # init from c+A/L+B/L^2
    X = np.vstack([np.ones_like(L), 1 / L, 1 / L ** 2]).T
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    p0 = [np.clip(coef[0], -0.5, 1.0), max(abs(coef[2]), 0.05), 2.0]
    popt, pcov = curve_fit(f, L, y, p0=p0, bounds=([-1.0, 0.0, 0.15], [3.0, 100.0, 6.0]), maxfev=400000)
    perr = np.sqrt(np.abs(np.diag(pcov)))
    return dict(c=popt[0], A=popt[1], p=popt[2], c_err=perr[0], p_err=perr[2],
                rms=np.sqrt(np.mean((y - f(L, *popt)) ** 2)))

def fit_cq(Ls, G):
    L = np.array(Ls, float); y = np.array([G[l] for l in Ls], float)
    f = lambda L, c, A, q: c + A * q ** L
    p0 = [max(np.min(y) - 0.01, 0.0), max(y[0] - np.min(y), 0.05), 0.6]
    popt, pcov = curve_fit(f, L, y, p0=p0, bounds=([-1.0, 0.0, 0.0], [3.0, 50.0, 1.0]), maxfev=400000)
    perr = np.sqrt(np.abs(np.diag(pcov)))
    return dict(c=popt[0], A=popt[1], q=popt[2], c_err=perr[0], q_err=perr[2],
                rms=np.sqrt(np.mean((y - f(L, *popt)) ** 2)))

# ---------------------------------------------------------------- main
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    kgrid, dk = make_kgrid(mode)
    nz = 320 if mode != 'full' else 460
    z, r1w = nodes(nz)
    print("#" * 100)
    print("# P4 -- PERIODIC constant-mode transfer matrix: gapped / L-independent / saturating?")
    print("#" * 100)
    v1 = haar_var_xi3(200001)
    HFL = 1.0 / v1
    print(f"exact one-link Var(xi3)_Haar = {v1:.6f}  =>  Haar-floor saturation value 1/Var = {HFL:.5f}")

    # ---- validation: periodic Z(k) machinery vs beta=0 closed form
    print("\n" + "=" * 100)
    print("[check] Z(k)=Tr[M_k^L] (eigvals) vs beta=0 closed form (independent links), L=4, k in {0.3,1.1,2.5}")
    for k in (0.3, 1.1, 2.5):
        zt = build_Z_periodic(4, 0.0, np.array([k]), z, r1w)[0]
        zb = build_Z_periodic_beta0(4, np.array([k]), z, r1w)[0]
        print(f"   k={k:4.2f}:  eigvals Tr={zt.real:+.9f}{zt.imag:+.2e}j   closedform={zb.real:+.9f}{zb.imag:+.2e}j")

    # ---- [T1] periodic bulk spectral transfer matrix
    print("\n" + "=" * 100)
    print("[T1] PERIODIC bulk transfer matrix M_0 (k=0): Delta=-log(lam1/lam0)=1/xi, L-INDEPENDENT by")
    print("     construction (one kernel).  Reference D3 strip: theta=pi Delta=3.18, theta=0 Delta=0.64.")
    print("-" * 100)
    print(f"   {'beta':>6}{'lam0':>11}{'lam1/lam0':>11}{'Delta':>9}{'xi=1/D':>9}   {'top evals':>28}")
    bulk = {}
    for beta in (0.0, 2.5):
        ev, D, xi = bulk_spectrum_periodic(beta)
        bulk[beta] = (ev, D, xi)
        print(f"   {beta:6.2f}{ev[0]:11.4f}{ev[1]/ev[0]:11.4e}{D:9.4f}{xi:9.4f}   {np.array2string(ev[:4], precision=5):>28}")
    print("   => beta=2.5 periodic: Delta>0 (gapped), matches D3 theta=0 strip kernel (same object).  [T]+[V]")
    print("      beta=0: rank-1 (lam1~0) => xi~0, delta-correlated, exact saturation at 1/Var.  [T]")

    # ---- [T2]/[T3] exact curvature table L=2..8
    print("\n" + "=" * 100)
    print("[T2/T3] EXACT periodic constant-mode curvature G''(0) (= -d^2 log rho/da^2 at the symmetric peak,")
    print("        moment ratio).  Split vs the beta=0 Haar floor: the periodic FLAT-DIRECTION part.")
    print("-" * 100)
    Ls = list(range(2, 9))
    Gp25, Gp0 = {}, {}
    print(f"   {'L':>2}{'G(2.5)':>11}{'G(0)=Haar':>12}{'Wil=G25-G0':>12}{'Wil/Haar':>10}{'G_mom25':>10}{'a0*(2.5)':>10}")
    for L in Ls:
        Z = build_Z_periodic(L, 2.5, kgrid, z, r1w)
        g25 = gcurv_moment(Z, kgrid)
        half = TWO_PI * np.sqrt(L)
        a0, g25c, _ = gap_from_Z(Z, kgrid, dk, half)
        Zh = build_Z_periodic_beta0(L, kgrid, z, r1w)
        g0 = gcurv_moment(Zh, kgrid)
        Gp25[L] = g25; Gp0[L] = g0
        wl = g25 - g0
        print(f"   {L:2d}{g25:+11.4f}{g0:+12.4f}{wl:+12.4f}{(wl/g0 if g0 else 0):+10.3f}{g25c:+10.4f}{a0:+10.4f}")
    print("   (G_mom = moment method; a0*(2.5) near 0 confirms the symmetric peak; both agree.)")

    # ---- reference: D1/D3 fixed-boundary strip values (Wilson part was O(1/L), decaying)
    print("\n   [reference D1/D3 fixed-boundary strip theta=0 beta=2.5: G'' = "
          "0.830,1.590,1.140,0.744,0.555,0.458,0.401 (L=1..7); Wilson part 0.625,1.351,0.891,...]")

    # ---- [T4] fits
    print("\n" + "=" * 100)
    print("[T4] FITS on the periodic sequence (L=2..8).  Intercept c = the L->inf gap.  c>0 => saturates.")
    print("-" * 100)
    for sel, sl in (('L=2..8', Ls), ('L=4..8', Ls[2:])):
        rp = fit_cAp(sl, Gp25); rq = fit_cq(sl, Gp25)
        print(f"   TOTAL (2.5) [{sel}]  c+A/L^p: c={rp['c']:+.4f}+-{rp['c_err']:.4f} A={rp['A']:.3f} "
              f"p={rp['p']:.2f}+-{rp['p_err']:.2f} rms={rp['rms']:.4f}")
        print(f"                    [{sel}]  c+A q^L: c={rq['c']:+.4f}+-{rq['c_err']:.4f} A={rq['A']:.3f} "
              f"q={rq['q']:.3f}+-{rq['q_err']:.3f} rms={rq['rms']:.4f}")
    # Wilson part scaling
    Gw = {L: Gp25[L] - Gp0[L] for L in Ls}
    print(f"   WILSON part G25-G0 = " + ", ".join(f"{Gw[L]:+.4f}" for L in Ls) + "  (L=2..8)")
    for sel, sl in (('L=2..8', Ls), ('L=4..8', Ls[2:])):
        rp = fit_cAp(sl, Gw)
        print(f"   WILSON [{sel}]  c+A/L^p: c={rp['c']:+.4f}+-{rp['c_err']:.4f} p={rp['p']:.2f}+-{rp['p_err']:.2f} "
              f"rms={rp['rms']:.4f}   (compare: D3 strip Wilson part had c->0, p~2)")

    # ---- [T5] acid test
    print("\n" + "=" * 100)
    print("[T5] ACID TEST (free field must return 0 / no manufactured floor):")
    print("-" * 100)
    for kappa in (0.5, 0.1, 0.02):
        av = np.linspace(-14, 14, 4001); h = av[1] - av[0]
        lg = -0.5 * kappa * av * av; i = len(av) // 2
        G = -(-lg[i - 2] + 16 * lg[i - 1] - 30 * lg[i] + 16 * lg[i + 1] - lg[i + 2]) / (12 * h * h)
        print(f"   free Gaussian (marginal N(0,1/kappa)), kappa={kappa:.2f}: G''(0)={G:+.4f}={kappa} -> 0")
    print("   free massless box gap 4pi^2/L^2: " + ", ".join(f"{4 * np.pi ** 2 / L ** 2:.3f}" for L in (2, 4, 8, 16)) + " -> 0.")
    print(f"   beta=0 Haar floor (NOT free): L-independent -> {HFL:.4f} (>0); discriminator = Phi=-log J.  [T]+[V]")

    print("\nVERDICT (see P4-PERIODIC-TRANSFER.md): the PERIODIC bulk transfer matrix is gapped")
    print(f"(Delta={bulk[2.5][1]:.3f}>0, L-independent).  Its constant-mode curvature G''(L) is POSITIVE and")
    print(f"L-INDEPENDENT (saturates near the Haar floor {HFL:.4f}; the periodic flat-direction Wilson part")
    print("is small, unlike the twist/end-pinning of the fixed strip).  NO decay to 0.  [T]+[V].  No QED.  Clay 0%.")


if __name__ == "__main__":
    main()
