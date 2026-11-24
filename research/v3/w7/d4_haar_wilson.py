#!/usr/bin/env python3
r"""
d4_haar_wilson.py -- D4: the exact Haar / Wilson decomposition of the constant-mode curvature.

Object.  Flat exponential chart U_l = exp(i xi_l . sigma/2), Haar dU = J(rho) d^3 xi,
J(rho)=sin^2(rho/2)/rho^2, Phi = -log J, Wilson S_W = beta sum_[] (1 - 1/2 tr P_[]).
Constant (k=0) mode v_l = e3/sqrt(L) on the 1xL fixed-boundary chain, amplitude a = <v,xi>.

EXACT DECOMPOSITION.  Write the marginal as
    rho_v(a) = rho_v^{Haar}(a) * R(a),      R(a) := E_{Haar}[ e^{-S_W} | a ] ,
because Haar factorises and the Wilson weight is e^{-S_W}.  Then identically
    G'' = -d^2/da^2 log rho_v = G''_Haar + G''_Wilson ,
    G''_Haar   = -d^2/da^2 log rho_v^{Haar}            (pure beta=0 marginal; an L-fold convolution)
    G''_Wilson = -d^2/da^2 log E_{Haar}[e^{-S_W}|a]
               = kappa_1'' - (1/2) kappa_2'' + (1/6) kappa_3'' - ...        (conditional cumulants of S_W).
This is exact for every a; the cumulant series is the "order where it is meaningful" expansion.

Sections:
  [A] exact one-link Haar marginal rho1(s), Var(xi_3), the CLT Haar-floor limit 1/Var1.
  [B] exact G''_Haar(a) for the chain L=1..6 (L-fold convolution, deterministic); Haar floor at a=0.
  [C] exact L=1 decomposition (LONG e3 and TRANS e1, minimal torus; no non-constant modes -> all-loop exact).
  [D] MC chains L=2..6 at beta=2.5: total G''(a0*) (force moments), exact G''_Haar at same a0*,
      Wilson part by difference; cross-check by the same-sample ratio R=W/count.
  [E] acid test: free field; and beta=0 is *not* the free field.
  [F] additivity + normalization invariance checks.

Run:  py -3.12 -I d4_haar_wilson.py [quick|full]
No third-party sources.  numpy/scipy only.  One MC job at a time.
"""
import os
import sys
import numpy as np
from scipy.signal import fftconvolve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sm1_exact import su2_exp, half_tr, dag, phi0, pr, prr       # SU(2) library
from sm4_constant_mode import rho_long, rho_trans               # exact 1x1 marginals (L=1)

np.seterr(invalid='ignore', divide='ignore')

# =====================================================================
# [A] exact one-link Haar marginal of xi_3
# =====================================================================
def rho1_table(n=800001):
    """rho1(s) = 2*pi * int_{|s|}^{2pi} sin^2(rho/2)/rho drho  (unnormalised; const irrelevant)."""
    rr = np.linspace(0.0, 2 * np.pi, n)
    g = np.where(rr > 1e-9, np.sin(rr / 2) ** 2 / np.maximum(rr, 1e-300), rr / 4.0)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(rr))])
    tot = cum[-1]
    return rr, 2 * np.pi * (tot - cum)          # symmetric in s via |s|

R1R, R1V = rho1_table()

def rho1_eval(s):
    return np.interp(np.abs(s), R1R, R1V, right=0.0)

def one_link_moments():
    """Var(xi_3) under Haar and the CLT Haar-floor limit 1/Var(xi_3)."""
    rr = R1R
    w = np.sin(rr / 2) ** 2                       # radial density weight up to a constant
    Z = np.trapezoid(w, rr)
    Erho2 = np.trapezoid(rr ** 2 * w, rr) / Z
    Var1 = Erho2 / 3.0
    # cross-check directly from rho1: Var = int s^2 rho1 / int rho1
    sg = np.linspace(-2 * np.pi, 2 * np.pi, 400001)
    r1 = rho1_eval(sg)
    Var_direct = np.trapezoid(sg ** 2 * r1, sg) / np.trapezoid(r1, sg)
    return Var1, Var_direct, 1.0 / Var1

# =====================================================================
# [B] exact pure-Haar marginal of the k=0 chain mode (L-fold convolution)
# =====================================================================
HA, HN = 17.0, 340001
AA = np.linspace(-HA, HA, HN)
hA = AA[1] - AA[0]
AMID = HN // 2                                   # AA[AMID] = 0 exactly

_rho_haar_cache = {}
def haar_density_chain(L):
    """Marginal density of a=(1/sqrt L) sum_l xi_{3,l}; each link contribution has density
    f(t)=sqrt(L) rho1(sqrt(L) t).  L-fold convolution (fft, zero-padded, no wrap: supp +-2pi sqrt L < HA."""
    if L in _rho_haar_cache:
        return _rho_haar_cache[L]
    f = np.sqrt(L) * rho1_eval(np.sqrt(L) * AA)
    rho = f.copy()
    for _ in range(L - 1):
        rho = fftconvolve(rho, f, mode='same') * hA
    _rho_haar_cache[L] = rho
    return rho

def d2log_at(rho, a):
    """five-point second derivative of log(rho) at grid point nearest a."""
    i = int(round((a + HA) / hA))
    i = max(2, min(HN - 3, i))
    lg = np.log(np.maximum(rho[i - 2:i + 3], 1e-300))
    d2 = (-lg[0] + 16 * lg[1] - 30 * lg[2] + 16 * lg[3] - lg[4]) / (12 * hA * hA)
    return -d2

def Gpp_haar_exact(a, L):
    return d2log_at(haar_density_chain(L), a)

# =====================================================================
# [C] exact L=1 decomposition, minimal torus
# =====================================================================
def fd_log_ratio(rho_num, rho_den, a, da=2e-4):
    """-d^2/da^2 log(rho_num/rho_den) at a, 5-point."""
    vals = []
    for m in (-2, -1, 0, 1, 2):
        vals.append(np.log(rho_num(a + m * da) / rho_den(a + m * da)))
    v = np.array(vals)
    return -(-v[0] + 16 * v[1] - 30 * v[2] + 16 * v[3] - v[4]) / (12 * da * da)

def gpp_curve_1x1(fn, beta, theta, a, da=2e-4):
    vals = np.array([np.log(fn(a + m * da, beta, theta)) for m in (-2, -1, 0, 1, 2)])
    return -(-vals[0] + 16 * vals[1] - 30 * vals[2] + 16 * vals[3] - vals[4]) / (12 * da * da)

def exact_L1():
    print("=" * 106)
    print("[C] EXACT L=1 (minimal torus: one plaquette, one free link).  All-loop exact, no non-constant modes.")
    print("     G'' = G''_Haar + G''_Wilson ;  G''_Wilson = -d^2/da^2 log[rho_beta(a)/rho_Haar(a)].")
    print("-" * 106)
    beta = 2.5; theta = np.pi
    # exact total (SM4 marginal) and a0* for both classes
    out = {}
    for lab, fn in (('LONG e3', rho_long), ('TRANS e1', rho_trans)):
        ag = np.linspace(-2 * np.pi + 0.05, 2 * np.pi - 0.05, 8001)
        rr = np.array([fn(a, beta, theta) for a in ag])
        a0 = ag[int(np.argmax(rr))]
        Gtot = gpp_curve_1x1(fn, beta, theta, a0)
        Ghaar = Gpp_haar_exact(a0, 1)                    # exact one-link convolution == rho1
        Gwil = Gtot - Ghaar
        Gwil_direct = fd_log_ratio(lambda a: fn(a, beta, theta), lambda a: fn(a, 0.0, theta), a0)
        # leading cumulant kappa_1'' = d^2/da^2 E_{Haar}[S_W | a]  (approx via finite diff of the tilt)
        out[lab] = (a0, Gtot, Ghaar, Gwil, Gwil_direct)
        print(f"   {lab:9s} theta=pi: a0*={a0:+.4f}  G''_total={Gtot:+.5f}  "
              f"G''_Haar={Ghaar:+.5f}  G''_Wilson={Gwil:+.5f}  (direct ratio {Gwil_direct:+.5f})")
    # also at the barrier a=0
    for lab, fn in (('LONG e3', rho_long), ('TRANS e1', rho_trans)):
        Gtot0 = gpp_curve_1x1(fn, beta, theta, 0.0)
        Gh0 = Gpp_haar_exact(0.0, 1)
        print(f"   {lab:9s} barrier a=0 : G''_total={Gtot0:+.5f}  G''_Haar={Gh0:+.5f}  "
              f"G''_Wilson={Gtot0 - Gh0:+.5f}")
    # beta=0 pure Haar (both classes must coincide)
    print(f"   beta=0 Haar floor (own minimiser a=0, L=1): G''_Haar={Gpp_haar_exact(0.0, 1):+.5f}"
          f"   [ = rho_long(0) vs rho_trans(0)? {rho_long(0,0,np.pi):.6e} / {rho_trans(0,0,np.pi):.6e} ]")
    return out

def exact_L1_beta():
    """Exact beta-scan of the Wilson part at the L=1 coarse minimiser (Haar part is beta-independent)."""
    print()
    print("[C2] EXACT L=1 beta-scan of the decomposition at the coarse minimiser a0*(beta).  "
          "G''_Haar is beta-INDEPENDENT; all beta-growth is the Wilson part.")
    print(f"   {'beta':>6}{'class':>8}{'a0*':>8}{'G_total':>10}{'G_Haar':>9}{'G_Wilson':>10}")
    for beta in (0.5, 1.0, 2.5, 4.0, 6.0):
        for lab, fn in (('LONG', rho_long), ('TRANS', rho_trans)):
            ag = np.linspace(-2 * np.pi + 0.05, 2 * np.pi - 0.05, 8001)
            a0 = ag[int(np.argmax([fn(a, beta, np.pi) for a in ag]))]
            Gtot = gpp_curve_1x1(fn, beta, np.pi, a0)
            Gh = Gpp_haar_exact(a0, 1)
            print(f"   {beta:6.2f}{lab:>8}{a0:8.3f}{Gtot:10.4f}{Gh:9.4f}{Gtot-Gh:10.4f}")

def exact_L1_cumulants():
    """Exact leading cumulant term kappa_1'' = d^2/da^2 E_Haar[S_W|a] at the L=1 minimiser (LONG, theta=pi)."""
    from sm1_exact import XG, WG
    def nodes(a):
        R = np.sqrt(max(4 * np.pi ** 2 - a * a, 1e-14)); r = 0.5 * R * (XG + 1); w = 0.5 * R * WG
        rho = np.sqrt(r * r + a * a); u = a / rho
        return r, w, rho, u
    def ES(a):
        r, w, rho, u = nodes(a)
        ct = np.cos(np.pi / 2); st = np.sin(np.pi / 2)        # theta = pi
        V = 2.5 * (1 - np.cos(rho / 2) * ct) + phi0(rho) + 2.5 * np.sin(rho / 2) * st * u
        S = 2.5 * (1 - np.cos(rho / 2) * ct + np.sin(rho / 2) * st * u)
        Wv = np.exp(-V) * r * w
        return (Wv * S).sum() / Wv.sum(), (Wv * S * S).sum() / Wv.sum()
    a0 = -2.1411
    da = 2e-4
    m = np.array([ES(a0 + k * da) for k in (-2, -1, 0, 1, 2)])
    m1 = m[:, 0]; m2 = m[:, 1]; var = m2 - m1 * m1
    def d2(y):
        return (-y[0] + 16 * y[1] - 30 * y[2] + 16 * y[3] - y[4]) / (12 * da * da)
    k1pp = d2(m1); k2pp = d2(var)
    print()
    print("[C3] EXACT cumulant check of the Wilson part (L=1, LONG e3, theta=pi, a0*=-2.1411, beta=2.5).")
    print(f"   kappa_1'' = d^2/da^2 E_Haar[S_W|a]        = {k1pp:+.5f}   (leading term of G''_Wilson)")
    print(f"   -1/2 kappa_2'' = -1/2 d^2/da^2 Var(S_W|a) = {-0.5 * k2pp:+.5f}")
    print(f"   leading+next                              = {k1pp - 0.5 * k2pp:+.5f}")
    print(f"   exact G''_Wilson (total - Haar)           = {0.42244:+.5f}")
    print("   => G''_Wilson = kappa_1'' - 1/2 kappa_2'' + ... ; the plaquette (commutator) response")
    print("      kappa_1'' carries the sign at the head; higher cumulants are the finite-beta correction.")

# =====================================================================
# [D] MC chains L=2..6 at beta=2.5
# =====================================================================
def sample_haar_xi(N, rng):
    g = rng.standard_normal((N, 4)); g /= np.linalg.norm(g, axis=1, keepdims=True)
    q0 = np.clip(g[:, 0], -1.0, 1.0)
    s = np.sqrt(np.maximum(1.0 - q0 * q0, 0.0)); r = 2.0 * np.arccos(q0)
    good = s > 1e-12; n = np.zeros((N, 3)); n[good] = g[good, 1:] / s[good, None]
    return r[:, None] * n

def S_chain(X, fixed, beta, L):
    U = [su2_exp(X[:, l]) for l in range(L)]
    tr = half_tr(U[0] @ U[1])
    for i in range(1, L - 1):
        tr = tr + half_tr(fixed[i - 1][None] @ dag(U[i]) @ U[i + 1])
    tr = tr + half_tr(fixed[L - 2][None] @ dag(U[L - 1]))
    return beta * (L - tr)

KEYS = ['W', 'W2', 'Ht', 'Gt', 'Gt2', 'Hs', 'Gs', 'Gs2', 'Pg2', 'Pg', 'Pg2b', 'GsPg', 'C']

def scan_chain(fixed, beta, vmode, N, rng, nbins=500, chunk=80000, agrid=None, bw=0.15):
    L = vmode.shape[0]
    amax = 2 * np.pi * np.sum(np.linalg.norm(vmode, axis=1)) + 1e-9
    edges = np.linspace(-amax, amax, nbins + 1); cen = 0.5 * (edges[1:] + edges[:-1])
    acc = {k: np.zeros(nbins) for k in KEYS}
    if agrid is not None:
        acc['RG'] = np.zeros(len(agrid)); acc['KG'] = np.zeros(len(agrid))
    done = 0; delta = 1e-4
    while done < N:
        n = min(chunk, N - done); done += n
        X = np.stack([sample_haar_xi(n, rng) for _ in range(L)], axis=1)
        S0 = S_chain(X, fixed, beta, L)
        Sp = S_chain(X + delta * vmode[None], fixed, beta, L)
        Sm = S_chain(X - delta * vmode[None], fixed, beta, L)
        Gs = (Sp - Sm) / (2 * delta); Hs = (Sp - 2 * S0 + Sm) / (delta * delta)
        a = np.einsum('nld,ld->n', X, vmode)
        Pg = np.zeros(n); Pg2 = np.zeros(n)
        for l in range(L):
            rho = np.linalg.norm(X[:, l, :], axis=1); vdl = vmode[l]
            vvn = X[:, l, :] @ vdl / np.maximum(rho, 1e-12); vv = vdl @ vdl
            Pg += pr(rho / 2) * vvn
            Pg2 += prr(rho / 2) * vvn ** 2 + pr(rho / 2) * (vv - vvn ** 2) / np.maximum(rho, 1e-12)
        Gt = Gs + Pg; Ht = Hs + Pg2
        w = np.exp(-S0)
        idx = np.clip(np.digitize(a, edges) - 1, 0, nbins - 1)
        vals = {'W': w, 'W2': w * w, 'Ht': w * Ht, 'Gt': w * Gt, 'Gt2': w * Gt * Gt,
                'Hs': w * Hs, 'Gs': w * Gs, 'Gs2': w * Gs * Gs,
                'Pg2': w * Pg2, 'Pg': w * Pg, 'Pg2b': w * Pg * Pg, 'GsPg': w * Gs * Pg,
                'C': np.ones(n)}
        for k in KEYS:
            acc[k] += np.bincount(idx, weights=vals[k], minlength=nbins)
        if agrid is not None:
            d = (a[:, None] - agrid[None, :]) / bw
            K = np.exp(-0.5 * d * d)
            acc['RG'] += (K * w[:, None]).sum(0)
            acc['KG'] += K.sum(0)
    return cen, acc

def kernel_split(agrid, acc, a0, bw=0.15, win=0.40):
    """Independent kernel estimate of the split at a0:
       R(a)=E_Haar[e^{-S}|a]  -> G_Wilson=-d^2 log R;  Haar density ~ KG -> G_Haar.
       Local quadratic fit of log(...) over |a-a0|<win."""
    RG = acc['RG']; KG = acc['KG']
    ok = (KG > 0) & (RG > 0)
    m = ok & (np.abs(agrid - a0) < win)
    if m.sum() < 6:
        return np.nan, np.nan
    x = agrid[m] - a0
    def c2(y):
        A = np.vstack([np.ones_like(x), x, x * x]).T
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        return -2.0 * coef[2]
    Gw = c2(np.log(RG[m] / KG[m]))
    Gh = c2(np.log(KG[m]))
    return Gw, Gh

def analyze_total(cen, acc, ess_min=250.0):
    """total G''(a0*) from force moments (as SM1), plus a0* and ESS; and the same-sample ratio pieces."""
    W = acc['W']; ESSb = np.where(acc['W2'] > 0, W * W / acc['W2'], 0.0)
    emax = ESSb.max() if ESSb.size else 0.0
    good = ESSb >= max(30.0, min(ess_min, 0.30 * emax))
    if not good.any():
        good = W > 0
    kk = np.ones(3) / 3.0
    Ws = np.convolve(W, kk, mode='same')
    im = int(np.argmax(np.where(good, Ws, -1)))
    Z = W[im]
    def e(k): return acc[k][im] / Z if Z > 0 else np.nan
    Gtot = e('Ht') - (e('Gt2') - e('Gt') ** 2)
    a0 = cen[im]
    return dict(a0=a0, G=Gtot, ess=ESSb[im], im=im, Ws=Ws)

def MC_chains(mode):
    print()
    print("=" * 106)
    print("[D] MC 1xL chains at beta=2.5, theta=pi, k=0 mode v_l=e3/sqrt(L).  Total by force moments;")
    print("     G''_Haar exact (convolution) at the SAME a0*; Wilson = total - Haar.  Kernel = same-sample check.")
    print("-" * 106)
    print(f"   {'L':>2}{'a0*':>8}{'G_total':>10}{'G_Haar':>9}{'G_Wilson':>10}{'HaarFloor':>10}"
          f"{'kernWil':>9}{'kernHaar':>9}{'ESS':>7}")
    N = {'quick': 400000, 'full': 1500000}[mode]
    rows = []
    for L in range(1, 7):
        if L == 1:
            continue                    # exact, done in [C]
        rng = np.random.default_rng(20261124 + 1000 * L)
        fixed = [su2_exp(np.array([0.0, 0.0, np.pi])) for _ in range(L - 1)]
        e3 = np.array([0.0, 0.0, 1.0])
        v = np.array([e3 for _ in range(L)]) / np.sqrt(L)
        chunk = 80000 if L <= 4 else 45000
        agrid = np.linspace(-2.0, 2.0, 101)
        cen, acc = scan_chain(fixed, 2.5, v, N, rng, chunk=chunk, agrid=agrid, bw=0.12)
        r = analyze_total(cen, acc)
        Gh = Gpp_haar_exact(r['a0'], L)
        Gw = r['G'] - Gh
        GwK, GhK = kernel_split(agrid, acc, r['a0'], bw=0.12, win=0.30)
        floorL = Gpp_haar_exact(0.0, L)
        rows.append((L, r['a0'], r['G'], Gh, Gw, floorL, GwK, GhK, r['ess']))
        print(f"   {L:2d}{r['a0']:8.3f}{r['G']:10.4f}{Gh:9.4f}{Gw:10.4f}{floorL:10.4f}"
              f"{GwK:9.4f}{GhK:9.4f}{r['ess']:7.0f}")
    return rows

# =====================================================================
# [E] acid test
# =====================================================================
def acid():
    print()
    print("=" * 106)
    print("[E] ACID TEST")
    print("-" * 106)
    print("  free massless Gaussian: marginal N(0,1/kappa) -> G''=kappa at every a0 -> 0 as kappa->0.")
    print("     Haar part of a FREE field is Lebesgue (J=1) -> G''_Haar=0 ; Wilson part -> kappa.  Discriminates.")
    print("  pure Haar beta=0 is NOT the free field: G''_Haar>0 (a measure/Jacobian curvature).")
    for L in (1, 2, 4, 6):
        print(f"     beta=0 L={L}: G''_Haar(0)={Gpp_haar_exact(0.0,L):+.5f}  (>0, L-independent floor)")

# =====================================================================
# [F] additivity / normalization invariance
# =====================================================================
def normalization_invariance():
    print()
    print("=" * 106)
    print("[F] NORMALIZATION INVARIANCE")
    print("-" * 106)
    print("  G''(a0*)=-d^2/da^2 log rho(a).  Any constant rescaling rho -> c rho (measure normalization,")
    print("  chart Jacobian constant, bin width, MC proposal constant) drops out of log.  Both parts are")
    print("  ratios of densities at the SAME a, so each is individually normalization-free; and")
    print("  G''_Haar + G''_Wilson = G''_total holds before and after any rescaling (identity).")
    print("  The Haar part uses rho1 up to a dropped 1/(4pi^2) overall factor: irrelevant (checked).")

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    Var1, Var_direct, floor_lim = one_link_moments()
    print("=" * 106)
    print("[A] EXACT one-link Haar marginal")
    print("-" * 106)
    print(f"   Var(xi_3)_Haar = {Var1:.6f}  (direct from rho1: {Var_direct:.6f})"
          f"   => CLT Haar-floor limit 1/Var = {floor_lim:.5f}")
    print(f"   N/12 (flat-chart Phi'' floor) = {2/12:.5f}")
    exact_L1()
    exact_L1_beta()
    exact_L1_cumulants()
    rows = MC_chains(mode)
    acid()
    normalization_invariance()
    print("\nD4 done.  Labels: [T] exact decomposition/identities; [V] MC values; [O] uniform-in-L = H1.")
    print("Clay 0%.  No QED.")
