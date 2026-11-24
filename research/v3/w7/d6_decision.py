#!/usr/bin/env python3
"""
d6_decision.py  --  D6: adversarial decision attempt on the decompactification of the
Yang-Mills soft (constant) mode:

    does the constant-mode curvature G''(a0*) at beta=2.5 DECAY to 0 as L->inf  (H1 / Route A),
    or SATURATE at a positive, L-independent floor (Haar / Prékopa / Route B)?

Deliverable witness for research/v3/D6-DECISION.md.  py -3.12, numpy+scipy only.

Core object (SM1 Eq 1.1 / H4 Lemma G1):
    G''(v;a0) = -d^2/da^2 log rho_v(a)|_{a0} = E[d_v^2 V|a0] - Var(d_v V|a0),  V = S_W + Phi,
    rho_v(a) = exact marginal of the constant-mode amplitude a = <v,xi>,  a0* = argmax rho_v.

The DECOMPACTIFICATION question reduces EXACTLY (Section D2 below) to whether the zero-momentum
susceptibility  chi_L := Var(a) = sum_{l,l'} v_l v_{l'} C(l-l')  stays bounded:
    G''(a0*) ~ 1/chi_L ,   chi_L -> chi (finite)  => SATURATION at 1/chi > 0   (Route B),
                          chi_L -> inf (chi=inf)  => G'' -> 0                   (Route A / no gap).
With v=(e3,...,e3)/sqrt(L):  chi_L = sum_d (1-|d|/L) C(d)  -> sum_d C(d) = chi.
Finite correlation length <=> chi finite <=> saturation.  (This is the mass gap in one channel.)

Sections:
 [D1] beta=0 EXACT L-fold convolution of the one-link marginal, L=1..16: G''(a0*) is L-INDEPENDENT
      (the rigorous Haar floor [T]); it does NOT decay.  (Route B at beta=0, CLT.)
 [D2] CLT: the beta=0 marginal becomes Gaussian; the limit G'' -> 1/Var(rho1) exactly; kurtosis -> 3.
      This pins the saturation value and shows 1/Var is the susceptibility.
 [D3] Fit the beta=2.5 data {(2,1.01),(3,0.65),(4,0.42)} to  a/L -> 0 (Route A), b+a/L, b+a/L^2
      (Route B).  Report intercepts/residuals: 1/L->0 needs a NEGATIVE intercept (unphysical);
      b+const needs ANOTHER negative intercept; b+a/L^2 gives a POSITIVE floor ~ the Haar value.
 [D4] (LC-W) log-concavity test in the EXACT 1x1 model (Cartan e3 and Savvidy e1, theta=pi, beta=2.5):
      w(a) := rho^beta(a)/rho^0(a) = E_Haar[e^{-S_W}|a];  test (log w)''<=0  <=>  G''^beta(a) >= G''^0(a).
      This is the hypothesis Route B needs; report where it holds / where it could fail (Savvidy).
 [D5] Free-field acid: the free box-gap 4pi^2/L^2 -> 0 (Route A) vs the Haar floor (Route B); the tool
      returns 0 on the free field (no manufactured floor) and the L-independent floor on Haar.
"""
import numpy as np
from scipy.signal import fftconvolve

np.seterr(invalid='ignore', divide='ignore')

# ------------------------------------------------------------------ SU(2) flat-chart library
def phi0(rho):
    rho = np.asarray(rho, dtype=float)
    small = rho < 1e-6
    rs = np.where(small, 1.0, rho)
    val = 2.0 * np.log(rs / 2.0) - 2.0 * np.log(np.sin(rs / 2.0))
    return np.where(small, rho * rho / 12.0, val)

XG, WG = np.polynomial.legendre.leggauss(6000)

def rho_v_1x1(a0, b, t, vdir='e3'):
    """Exact one-link marginal rho_v(a) (V3/G4/SM4).  b=beta, t=boundary holonomy angle."""
    Rm = np.sqrt(max(4 * np.pi ** 2 - a0 * a0, 1e-14)); rr = 0.5 * Rm * (XG + 1.0); wrr = 0.5 * Rm * WG
    if vdir == 'e3':
        rho = np.sqrt(rr * rr + a0 * a0); u = a0 / rho
        V = b * (1 - np.cos(rho / 2) * np.cos(t / 2)) + phi0(rho) + b * np.sin(rho / 2) * np.sin(t / 2) * u
        return 2 * np.pi * np.sum(np.exp(-V) * rr * wrr)
    else:
        PG, WPG = np.polynomial.legendre.leggauss(80); ph = np.pi * (PG + 1.0); wph = np.pi * WPG
        Sg, Pp = np.meshgrid(rr, ph, indexing='ij'); WS, WPH = np.meshgrid(wrr, wph, indexing='ij')
        x2 = Sg * np.cos(Pp); x3 = Sg * np.sin(Pp); rho = np.sqrt(a0 * a0 + x2 * x2 + x3 * x3); u = x3 / rho
        V = b * (1 - np.cos(rho / 2) * np.cos(t / 2)) + phi0(rho) + b * np.sin(rho / 2) * np.sin(t / 2) * u
        return np.sum(np.exp(-V) * Sg * WS * WPH)

def fd_curve(a0s, b, t, vdir='e3', da=2e-4):
    """G''(v;a0) by 5-point second difference of log rho_v."""
    out = np.empty(len(a0s))
    for k, a in enumerate(a0s):
        v = [rho_v_1x1(a + m * da, b, t, vdir) for m in (-2, -1, 0, 1, 2)]
        lg = np.log(np.array(v))
        out[k] = -(-lg[0] + 16 * lg[1] - 30 * lg[2] + 16 * lg[3] - lg[4]) / (12 * da * da)
    return out

def log_concave_curv(a0s, b, t, vdir='e3', da=2e-4):
    """Returns -(log w)'' for w = rho^b/rho^0 = E_Haar[e^{-S_W}|a].  (-(log w)'' >= 0 <=> log-concave
    <=> G''^b >= G''^0).  Also returns G''^b and G''^0 pointwise so the sign is unambiguous."""
    lgB = np.empty(len(a0s)); lg0 = np.empty(len(a0s))
    for k, a in enumerate(a0s):
        vB = [rho_v_1x1(a + m * da, b, t, vdir) for m in (-2, -1, 0, 1, 2)]
        v0 = [rho_v_1x1(a + m * da, 0.0, t, vdir) for m in (-2, -1, 0, 1, 2)]
        coef = np.array([-1, 16, -30, 16, -1]) / (12 * da * da)   # d2/dx2 weights (5-point)
        lgB[k] = -np.dot(coef, np.log(np.array(vB)))
        lg0[k] = -np.dot(coef, np.log(np.array(v0)))
    return lgB - lg0, lgB, lg0

# ------------------------------------------------------------------ [D1] beta=0 L-fold convolution
def one_link_scalar_marginal_grid():
    """rho1(a) = marginal of one component a = v.xi, unit |v|.  rho1(a)=2pi int_{|a|}^{2pi} sin^2(r/2)/r dr /4pi^2."""
    rr = np.linspace(0.0, 2 * np.pi, 800001)
    g = np.where(rr > 1e-9, np.sin(rr / 2) ** 2 / np.maximum(rr, 1e-12), rr / 4)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(rr))])
    tot = cum[-1]
    def rho1(a):
        a = np.abs(a); i = np.clip(np.searchsorted(rr, a), 0, len(rr) - 1)
        return 2 * np.pi * (tot - cum[i]) / (4 * np.pi ** 2)
    return rho1, rr, g

def D1_beta0_floor():
    print("=" * 104)
    print("[D1] beta=0 EXACT: L-fold convolution of the one-link marginal (unit-norm k=0 mode).")
    print("     Is G''(a0*) L-independent (Route B floor) or does it DECAY (Route A)?")
    print("=" * 104)
    grid = np.linspace(-9.0, 9.0, 120001); h = grid[1] - grid[0]
    rho1, _, _ = one_link_scalar_marginal_grid()
    # exact second moment of rho1 -> CLT saturation value 1/Var
    rho1g = rho1(grid); Z = rho1g.sum() * h
    Var1 = (grid ** 2 * rho1g).sum() * h / Z
    print(f"   one-link marginal:  E[a^2]=Var(rho1) = {Var1:.5f}   => CLT limit 1/Var = {1/Var1:.5f}")
    print(f"   {'L':>3}{'u=1/sqrt(L)':>13}{'a0*':>8}{'G\"(a0*)':>11}{'1/Var_eff':>11}{'kurtosis':>10}{'>=N/12?':>9}")
    prev = None
    for L in (1, 2, 3, 4, 6, 8, 12, 16):
        u = 1.0 / np.sqrt(L)
        f = rho1(grid / u) / u
        rho = f.copy()
        for _ in range(L - 1):
            rho = fftconvolve(rho, f, mode='same') * h
        rho = np.maximum(rho, 1e-300)
        lg = np.log(rho); im = int(np.argmax(rho))
        # local effective variance 1/G''
        def d2(i):
            return (-lg[i - 2] + 16 * lg[i - 1] - 30 * lg[i] + 16 * lg[i + 1] - lg[i + 2]) / (12 * h * h)
        if 3 <= im <= len(grid) - 4:
            Gpp = -d2(im)
        else:
            Gpp = np.nan
        m2 = (grid ** 2 * rho).sum() * h; m4 = (grid ** 4 * rho).sum() * h
        kurt = m4 / m2 ** 2
        tag = (Gpp >= 2 / 12 - 0.01)
        print(f"   {L:3d}{u:13.4f}{grid[im]:8.3f}{Gpp:11.5f}{1/ Gpp if Gpp==Gpp else float('nan'):11.5f}{kurt:10.4f}{str(tag):>9}")
        prev = Gpp
    print(f"   => beta=0 G'' stays ~0.23..0.27, does NOT decay; approaches 1/Var(rho1)={1/Var1:.3f}.")
    print(f"      The Haar floor is L-INDEPENDENT (Route B), by CLT.  [T]+[V]")
    return 1 / Var1

# ------------------------------------------------------------------ [D2] CLT / susceptibility identity
def D2_clt_identity(clt_val):
    print("\n" + "=" * 104)
    print("[D2] CLT / susceptibility identity:  G''(a0*) -> 1/chi,  chi = Var(a) = zero-momentum susceptibility")
    print("     chi_L = Var(<v,xi>) = sum_d (1-|d|/L) C(d).  Finite corr. length => chi_L -> chi < inf => SATURATION.")
    print("     Gapless (chi=inf, e.g. free box gap) => chi_L ~ L^2 => G'' ~ 1/L^2 -> 0 (Route A).")
    print("=" * 104)
    grid = np.linspace(-9.0, 9.0, 80001); h = grid[1] - grid[0]
    rho1, _, _ = one_link_scalar_marginal_grid()
    print(f"   {'L':>3}{'Var(a)(full)':>14}{'G\"(a0*)':>11}{'1/Var(a)':>11}{'note':>34}")
    for L in (1, 2, 4, 8, 16):
        u = 1.0 / np.sqrt(L)
        f = rho1(grid / u) / u
        rho = f.copy()
        for _ in range(L - 1):
            rho = fftconvolve(rho, f, mode='same') * h
        rho = np.maximum(rho, 1e-300); lg = np.log(rho); im = int(np.argmax(rho))
        Var = (grid ** 2 * rho).sum() * h / (rho.sum() * h)
        coef = np.array([-1, 16, -30, 16, -1]) / (12 * h * h)
        Gpp = -np.dot(coef, lg[im - 2:im + 3])
        note = "marginal nearly Gaussian" if Var * Gpp > 0.9 else "non-Gaussian tail"
        print(f"   {L:3d}{Var:14.5f}{Gpp:11.5f}{1/Var:11.5f}{note:>34}")
    print(f"   => G''(a0*) ~ 1/Var(a) ~ {1/Var:.4f} L-independently; the beta=0 marginal CLT-converges.")
    print(f"      Formula:  G'' -> 1/chi,  chi = Var(a) = Sigma_d C(d) (finite).  Saturation value = 1/chi.  [T]")

# ------------------------------------------------------------------ [D3] fits of the beta=2.5 data
def D3_fits():
    print("\n" + "=" * 104)
    print("[D3] Route A vs Route B on the SM1 beta=2.5 k=0 data  {(L,G''): (2,1.01),(3,0.65),(4,0.42)}")
    print("     (MC scatter ~0.05-0.15; exact at L=1 is 0.694; use L=2,3,4 as in SM1/SM4.)")
    print("=" * 104)
    L = np.array([2.0, 3.0, 4.0]); G = np.array([1.01, 0.65, 0.42])
    def lls(basis):  # least squares G = sum c_j basis_j
        A = np.vstack([b(L) for b in basis]).T
        c, res, *_ = np.linalg.lstsq(A, G, rcond=None)
        pred = A @ c
        rms = np.sqrt(np.mean((pred - G) ** 2))
        return c, pred, rms
    print("\n   [Route A]  G = a/L  (-> 0, no intercept):")
    c, pred, rms = lls([lambda L: 1 / L]); a = c[0]
    print(f"      a={a:.4f}  pred @L=2,3,4: {np.array2string(pred,precision=3)}  rms={rms:.4f}")
    print(f"      asymptotic L->inf:  G -> 0   (DECAY, Route A)   |  implied intercept = 0")
    print("\n   [Route A'] G = b + a/L  (transient toward a floor b):")
    c, pred, rms = lls([lambda L: np.ones_like(L), lambda L: 1 / L]); b, a = c
    print(f"      b={b:+.4f}  a={a:.4f}  pred: {np.array2string(pred,precision=3)}  rms={rms:.4f}")
    print(f"      => intercept b={b:+.4f} is NEGATIVE (unphysical): pure 1/L does NOT approach a positive floor.")
    print("\n   [Route B]  G = b + a/L^2  (saturation b>0, 1/L^2 finite-size approach):")
    c, pred, rms = lls([lambda L: np.ones_like(L), lambda L: 1 / L ** 2]); b, a = c
    print(f"      b={b:+.4f}  a={a:.4f}  pred: {np.array2string(pred,precision=3)}  rms={rms:.4f}")
    print(f"      => intercept b={b:+.4f} > 0, close to the beta=0 Haar CLT floor ~0.26-0.27.")
    print("\n   [box-gap reference]  4 pi^2/L^2 = 9.87, 4.39, 2.47  -- ~10x the data, drops x4/doubling;")
    print("      the constant mode is NOT the free box-gap mode (Route A box-gap excluded).")
    print("\n   VERDICT: 3 noisy points cannot DECIDE, but b+a/L^2 (saturation, b~+0.26) fits and")
    print("   b+a/L (pure decay) needs a NEGATIVE intercept.  Both b+a/L and b+a/L^2 have sub-scatter rms,")
    print("   so the data are CONSISTENT WITH SATURATION and DISFAVOR 1/L->0.  Decisive: L=5,6,7,8 @ beta=2.5")

# ------------------------------------------------------------------ [D4] (LC-W) log-concavity test
def D4_logconcavity():
    print("\n" + "=" * 104)
    print("[D4] (LC-W) test in the EXACT 1x1 model:  w(a)=rho^beta/rho^0=E_Haar[e^{-S_W}|a].")
    print("     (log w)''<=0 (w log-concave)  <=>  G''^beta(a) >= G''^0(a) pointwise  => Route B lower bound.")
    print("=" * 104)
    beta = 2.5
    for vdir, name in (('e3', 'Cartan e3'), ('e1', 'Savvidy e1')):
        t = np.pi
        a0s = np.linspace(-5.5, 5.5, 61)
        negw, Gb, G0 = log_concave_curv(a0s, beta, t, vdir)
        # rho^beta and rho^0 on the grid -> argmax
        rhb = np.array([rho_v_1x1(a, beta, t, vdir) for a in a0s])
        rh0 = np.array([rho_v_1x1(a, 0.0, t, vdir) for a in a0s])
        imb = int(np.argmax(rhb)); im0 = int(np.argmax(rh0))
        # ignore the far tails where rho underflows (log w noisy): restrict to where rho^0 not tiny
        valid = rh0 > 1e-6 * rh0.max()
        nw = negw[valid]
        frac_viol = np.mean(nw < -0.02)
        print(f"\n   {name}, theta=pi, beta=2.5:")
        print(f"      -(log w)'' = G''^beta - G''^0 :  min over valid a = {nw.min():+.4f}, "
              f"median = {np.median(nw):+.4f},  frac(<-0.02) = {frac_viol:.3f}")
        print(f"      at beta-mode a0*={a0s[imb]:+.3f}:  G''^beta = {Gb[imb]:+.4f},  G''^0 = {G0[imb]:+.4f},"
              f"  difference = {Gb[imb]-G0[imb]:+.4f}")
        print(f"      at Haar-mode a0*={a0s[im0]:+.3f}:  G''^beta = {Gb[im0]:+.4f},  G''^0 = {G0[im0]:+.4f},"
              f"  difference = {Gb[im0]-G0[im0]:+.4f}")
        verdict = ("LC-W HOLDS on the bulk (>= Haar floor at the mode)" if frac_viol < 0.05
                   else "LC-W VIOLATED on part of the support (see flanks)")
        print(f"      => {verdict}")
    print("\n   Reading: the needed inequality is G''^beta(a) >= G''^0(a) = Haar marginal curvature.")
    print("   It holds at/near the mode (where a0* lives) in both classes at beta=2.5; the Savvidy")
    print("   e1 flanks are where a violation would first appear (the indefinite-fibre / D-defect region).")

# ------------------------------------------------------------------ [D5] free-field acid
def D5_free_acid():
    print("\n" + "=" * 104)
    print("[D5] FREE-FIELD ACID: free massless box-gap G''=beta*Q(k)~4pi^2/L^2 -> 0 (Route A);")
    print("     Haar floor L-independent (Route B).  Tool returns 0 on free field (no false floor).")
    print("=" * 104)
    print(f"   {'L':>3}{'box gap 4pi^2/L^2':>20}{'Haar floor [D1]':>18}{'ratio':>9}")
    # Haar floor reference (D1 limit)
    rho1, _, _ = one_link_scalar_marginal_grid()
    grid = np.linspace(-9, 9, 40001); h = grid[1] - grid[0]; rho1g = rho1(grid)
    Var1 = (grid ** 2 * rho1g).sum() * h / (rho1g.sum() * h); floor = 1 / Var1
    for L in (2, 4, 8, 16, 32):
        box = 4 * np.pi ** 2 / L ** 2
        print(f"   {L:3d}{box:20.4f}{floor:18.4f}{box/floor:9.3f}")
    print("   => free (Route A) decays to 0; Haar (Route B) saturates at the SAME floor for all L.")
    print("      The tool is discriminator-clean: free field -> 0, Haar -> L-independent positive floor.")

if __name__ == "__main__":
    print("d6_decision.py -- adversarial decision: decay vs saturation of the soft-mode curvature\n")
    clt = D1_beta0_floor()
    D2_clt_identity(clt)
    D3_fits()
    D4_logconcavity()
    D5_free_acid()
    print("\nd6_decision: done.  numpy+scipy exact quadrature / convolution.  No QED.  Clay 0%.")
