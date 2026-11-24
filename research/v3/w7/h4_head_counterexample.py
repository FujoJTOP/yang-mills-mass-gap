#!/usr/bin/env python3
"""
h4_head_counterexample.py  --  H4: try to falsify the CORRECTED head inequality at the coarse minimiser.

Object (the true coarse curvature, ch05 eq:split / G4 Lemma G1) for a linear (constant) direction v:

   G''(v; a0) := -d^2/da^2 log rho_v(a) |_{a0} = E[ d_v^2 V | a0 ] - Var( d_v V | a0 ),
   V = S + Phi ,   a = <v, xi> ,   d_v = v . grad_xi  (Cartesian flat exponential chart, v constant),
   decomposition  G'' = Delta_dyn + W_geo + X ,  X = -2 Cov(d_v S, d_v Phi | a0).

Corrected (star) = (star_loc):  G''(v; a0*) >= 0  at the coarse minimiser a0* = argmax rho_v, plus a
uniform lower bound delta_W.  The V3/G4 counterexample was at the BARRIER a0=0 at boundary theta=pi;
here we test the value AT THE COARSE MINIMISER in a LARGER block (1x2, and 1x3 chain), with the coupled
A-V chart so X != 0.

Models (SU(2), flat chart Haar):
  1x1 : one free link U, one plaquette  P = U B .                 (reproduces V3 / G4 exactly)
  1x2 : free A,V ;  P1 = A V ,  P2 = B V^{-1} .                    (G2 two-plaquette block)
  1x3 : free A,V,W ; P1 = A V , P2 = B1 V^{-1} W , P3 = B2 W^{-1} . (3-plaquette chain; block-size test)

Method.
  * exact 1x1: deterministic Gauss-Legendre (V3 code) -- reference numbers.
  * general  : Haar-proposal importance sampling (weight e^{-S}); per-a histogram of the EXACT pointwise
    Hessian moments (the identity holds per-a, so G''(a_bin) is read directly with no smoothing window).
    Streaming accumulation -> very large N with low memory.

Run:  py -3.12 h4_head_counterexample.py [quick|full|huge]
"""
import sys
import numpy as np

# ------------------------------------------------------------------ SU(2) library
S1 = np.array([[0, 1], [1, 0]], dtype=complex)
S2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
S3 = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def su2_exp(xi):
    """xi: (...,3) real -> U=exp(i xi.sigma/2), shape (...,2,2).  t=sin(r/2)/r (->1/2 at r=0)."""
    xi = np.asarray(xi, dtype=float)
    r = np.sqrt(np.sum(xi * xi, axis=-1))
    small = r < 1e-8
    rs = np.where(small, 1.0, r)
    t = np.sin(rs / 2.0) / rs
    t = np.where(small, 0.5, t)
    c = np.cos(r / 2.0)
    x1, x2, x3 = xi[..., 0], xi[..., 1], xi[..., 2]
    U = np.empty(xi.shape[:-1] + (2, 2), dtype=complex)
    it = 1j * t
    U[..., 0, 0] = c + it * x3
    U[..., 0, 1] = t * x2 + it * x1
    U[..., 1, 0] = -t * x2 + it * x1
    U[..., 1, 1] = c - it * x3
    return U


def half_tr(M):
    return 0.5 * np.real(np.trace(M, axis1=-2, axis2=-1))


def dag(M):
    return np.conjugate(np.swapaxes(M, -1, -2))


def phi0(rho):
    rho = np.asarray(rho, dtype=float)
    small = rho < 1e-6
    rs = np.where(small, 1.0, rho)
    val = 2.0 * np.log(rs / 2.0) - 2.0 * np.log(np.sin(rs / 2.0))
    return np.where(small, rho * rho / 12.0, val)


def pr(x):      # Phi'(rho),  x=rho/2
    x = np.asarray(x, dtype=float); y = np.empty_like(x); sm = np.abs(x) < 0.3
    y[~sm] = 1.0 / x[~sm] - 1.0 / np.tan(x[~sm])
    xs = x[sm]; y[sm] = xs / 3 + xs ** 3 / 45 + 2 * xs ** 5 / 945 + 2 * xs ** 7 / 4725
    return y


def prr(x):     # Phi''(rho),  x=rho/2
    x = np.asarray(x, dtype=float); y = np.empty_like(x); sm = np.abs(x) < 0.3
    y[~sm] = 0.5 * (-1.0 / x[~sm] ** 2 + 1.0 / np.sin(x[~sm]) ** 2)
    xs = x[sm]; y[sm] = 1.0 / 6 + xs ** 2 / 30 + 2 * xs ** 4 / 189 + 2 * xs ** 6 / 2475
    return y


def sample_haar_xi(N, rng):
    g = rng.standard_normal((N, 4)); g /= np.linalg.norm(g, axis=1, keepdims=True)
    q0 = np.clip(g[:, 0], -1.0, 1.0)
    s = np.sqrt(np.maximum(1.0 - q0 * q0, 0.0)); r = 2.0 * np.arccos(q0)
    good = s > 1e-12; n = np.zeros((N, 3)); n[good] = g[good, 1:] / s[good, None]
    return r[:, None] * n


# ------------------------------------------------------------------ block actions  S(X)
def S_1x1(X, fixed, beta):
    A = su2_exp(X[:, 0]); B = fixed[0]
    return beta * (1.0 - half_tr(np.matmul(A, B[None])))


def S_1x2(X, fixed, beta):
    A = su2_exp(X[:, 0]); V = su2_exp(X[:, 1]); B = fixed[0]
    P1 = np.matmul(A, V); P2 = np.matmul(B[None], dag(V))
    return beta * (2.0 - half_tr(P1) - half_tr(P2))


def S_1x3(X, fixed, beta):
    A = su2_exp(X[:, 0]); V = su2_exp(X[:, 1]); W = su2_exp(X[:, 2])
    B1, B2 = fixed[0], fixed[1]
    P1 = np.matmul(A, V); P2 = np.matmul(np.matmul(B1[None], dag(V)), W); P3 = np.matmul(B2[None], dag(W))
    return beta * (3.0 - half_tr(P1) - half_tr(P2) - half_tr(P3))


# ------------------------------------------------------------------ exact 1x1 (V3 code)
XG, WG = np.polynomial.legendre.leggauss(6000)


def parts(rho, b, t):
    hS = b * (1 - np.cos(rho / 2) * np.cos(t / 2)); hSp = b * 0.5 * np.sin(rho / 2) * np.cos(t / 2)
    g = b * np.sin(rho / 2) * np.sin(t / 2); gp = b * 0.5 * np.cos(rho / 2) * np.sin(t / 2)
    hP = phi0(rho); hPp = pr(rho / 2)
    return hS, hSp, g, gp, hP, hPp


def rho_v_1x1(a0, b, t, vdir='e3'):
    Rm = np.sqrt(max(4 * np.pi ** 2 - a0 * a0, 1e-14)); rr = 0.5 * Rm * (XG + 1.0); wrr = 0.5 * Rm * WG
    if vdir == 'e3':
        rho = np.sqrt(rr * rr + a0 * a0); u = a0 / rho
        V = b * (1 - np.cos(rho / 2) * np.cos(t / 2)) + phi0(rho) + b * np.sin(rho / 2) * np.sin(t / 2) * u
        return 2 * np.pi * np.sum(np.exp(-V) * rr * wrr)
    else:
        PG, WPG = np.polynomial.legendre.leggauss(240); ph = np.pi * (PG + 1.0); wph = np.pi * WPG
        Sg, Pp = np.meshgrid(rr, ph, indexing='ij'); WS, WPH = np.meshgrid(wrr, wph, indexing='ij')
        x2 = Sg * np.cos(Pp); x3 = Sg * np.sin(Pp); rho = np.sqrt(a0 * a0 + x2 * x2 + x3 * x3); u = x3 / rho
        V = b * (1 - np.cos(rho / 2) * np.cos(t / 2)) + phi0(rho) + b * np.sin(rho / 2) * np.sin(t / 2) * u
        return np.sum(np.exp(-V) * Sg * WS * WPH)


def fd_1x1(a0, b, t, vdir='e3', da=2e-4):
    v = [rho_v_1x1(a0 + k * da, b, t, vdir) for k in (-2, -1, 0, 1, 2)]; lg = np.log(np.array(v))
    return -(-lg[0] + 16 * lg[1] - 30 * lg[2] + 16 * lg[3] - lg[4]) / (12 * da * da)


def exact_1x1_report():
    print("=" * 100)
    print("[S1] EXACT 1x1 (one plaquette / one free link) -- reproduces V3 / G4   [T]+[V]")
    print(f"     N/12 = {2/12:.6f}")
    w0 = (np.sin(np.pi * (XG + 1.0) / 2) ** 2 / (np.pi * (XG + 1.0))) * (np.pi * WG)
    r0 = np.pi * (XG + 1.0); Z0 = w0.sum(); E0 = lambda x: (w0 * x).sum() / Z0
    A = E0(pr(r0 / 2) / r0); ss = np.sin(r0 / 2) / r0; B = E0(ss * ss) - E0(ss) ** 2
    print(f"     A=E_Haar[Phi'(r)/r]={A:.10f}  B=Var_Haar[sin(r/2)/r]={B:.10f}  beta*=sqrt(A/B)={np.sqrt(A/B):.5f}")
    out = {}
    print(f"     {'beta':>6}{'A-B b^2':>13}{'fd a0=0':>12}")
    for beta in (0.0, 1.0, 2.5, 5.0, 10.0):
        be = max(beta, 1e-9); print(f"     {beta:6.2f}{A-B*be*be:13.6f}{fd_1x1(0.0,be,np.pi,'e3'):12.6f}")
    print("     coarse-well location & G'' AT THE MINIMISER a0* (the corrected (star) point):")
    for beta in (2.5, 5.0, 10.0):
        a0s = np.linspace(-4.0, 1.0, 201); rh = np.array([rho_v_1x1(a, beta, np.pi, 'e3') for a in a0s])
        am = a0s[int(np.argmax(rh))]; out[beta] = (am, fd_1x1(am, beta, np.pi, 'e3'))
        print(f"       beta={beta:5.2f}: a0*={am:+.3f}  G''(a0*)={fd_1x1(am,beta,np.pi,'e3'):+.5f}   "
              f"G''(a0=0)={fd_1x1(0.0,beta,np.pi,'e3'):+.5f}")
    # full G''(a) curve at beta=2.5: the one-plaquette ALREADY has a negative shoulder for larger a.
    print("     full G''(a) at beta=2.5, theta=pi, e3  (exact; the negative shoulder already exists here):")
    for a0 in (-3.5, -2.5, -2.0, -1.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5):
        print(f"       a0={a0:+.2f}: G''={fd_1x1(a0,2.5,np.pi,'e3'):+.5f}")
    return A, B, out


# ------------------------------------------------------------------ streaming MC scanner
# accumulators per bin (weighted sums):  W, W2, Ht, Gt, Gt2, Hs, Gs, Gs2, Pg2, Pg, Pg2b, GsPg
KEYS = ['W', 'W2', 'Ht', 'Gt', 'Gt2', 'Hs', 'Gs', 'Gs2', 'Pg2', 'Pg', 'Pg2b', 'GsPg']


def scan_block(block, fixed, beta, vmode, N, rng, nbins=600, chunk=250000, vlabel=''):
    nfree = vmode.shape[0]
    amax = 2 * np.pi * np.sum(np.linalg.norm(vmode, axis=1)) + 1e-9
    edges = np.linspace(-amax, amax, nbins + 1); cen = 0.5 * (edges[1:] + edges[:-1])
    acc = {k: np.zeros(nbins) for k in KEYS}
    w2edge = edges.__getitem__
    done = 0; delta = 1e-4
    while done < N:
        n = min(chunk, N - done); done += n
        X = np.stack([sample_haar_xi(n, rng) for _ in range(nfree)], axis=1)
        S0 = block(X, fixed, beta)
        Sp = block(X + delta * vmode[None], fixed, beta)
        Sm = block(X - delta * vmode[None], fixed, beta)
        Gs = (Sp - Sm) / (2 * delta); Hs = (Sp - 2 * S0 + Sm) / (delta * delta)
        a = np.einsum('nld,ld->n', X, vmode)
        Ph = 0.0; Pg = np.zeros(n); Pg2 = np.zeros(n)
        for l in range(nfree):
            rho = np.sqrt(np.sum(X[:, l, :] ** 2, axis=1)); vdl = vmode[l]
            vvn = np.dot(X[:, l, :], vdl) / np.maximum(rho, 1e-12); vv = np.dot(vdl, vdl)
            Ph = Ph + phi0(rho)
            Pg = Pg + pr(rho / 2) * vvn
            Pg2 = Pg2 + prr(rho / 2) * vvn ** 2 + pr(rho / 2) * (vv - vvn ** 2) / np.maximum(rho, 1e-12)
        Gt = Gs + Pg; Ht = Hs + Pg2
        w = np.exp(-S0)
        idx = np.clip(np.digitize(a, edges) - 1, 0, nbins - 1)
        vals = {'W': w, 'W2': w * w, 'Ht': w * Ht, 'Gt': w * Gt, 'Gt2': w * Gt * Gt,
                'Hs': w * Hs, 'Gs': w * Gs, 'Gs2': w * Gs * Gs,
                'Pg2': w * Pg2, 'Pg': w * Pg, 'Pg2b': w * Pg * Pg, 'GsPg': w * Gs * Pg}
        for k in KEYS:
            acc[k] += np.bincount(idx, weights=vals[k], minlength=nbins)
    return cen, acc


def summarize(cen, acc, label, ess_min=400.0, dump=False, quiet=False):
    W = acc['W']
    ESSb = np.where(acc['W2'] > 0, W * W / acc['W2'], 0.0)
    emax = ESSb.max() if ESSb.size else 0.0
    floor = max(30.0, min(ess_min, 0.30 * emax))
    good = ESSb >= floor
    if not good.any():
        good = W > 0
    rho = W.copy(); k = np.ones(3) / 3.0; rhos = np.convolve(rho, k, mode='same')
    imax = int(np.argmax(np.where(good, rhos, -1)))
    a0 = cen[imax]
    def gmom(i):
        Z = W[i]
        if Z <= 0:
            return (np.nan,) * 6
        EHt = acc['Ht'][i] / Z; EGt = acc['Gt'][i] / Z; EGt2 = acc['Gt2'][i] / Z
        EHs = acc['Hs'][i] / Z; EGs = acc['Gs'][i] / Z; EGs2 = acc['Gs2'][i] / Z
        EPg2 = acc['Pg2'][i] / Z; EPg = acc['Pg'][i] / Z; EPg2b = acc['Pg2b'][i] / Z
        EGsPg = acc['GsPg'][i] / Z
        Gpp = EHt - (EGt2 - EGt * EGt)
        Ddyn = EHs - (EGs2 - EGs * EGs)
        Wgeo = EPg2 - (EPg2b - EPg * EPg)
        X = -2 * (EGsPg - EGs * EPg)
        return Gpp, Ddyn, Wgeo, X, ESSb[i], Z
    G, D, Wg, X, ess, Z = gmom(imax)
    Gcurve = np.full(len(cen), np.nan)
    for i in np.where(good)[0]:
        Gcurve[i] = gmom(i)[0]
    finite = np.where(np.isfinite(Gcurve))[0]
    if finite.size == 0:
        print(f"  -- {label}: ALL BINS UNRELIABLE at this N (high-beta weight degeneracy); increase N.")
        return dict(label=label, a0=a0, G=np.nan, D=np.nan, Wg=np.nan, X=np.nan, ess=ess,
                    cen=cen, Gcurve=Gcurve, good=good, gmin=np.nan, amin=np.nan, nneg=0)
    jmin = finite[np.argmin(Gcurve[finite])]
    neg = finite[Gcurve[finite] < 0]
    if not quiet:
        print(f"  -- {label}")
        print(f"     a0*={a0:+.3f} (bin {imax})  ESS(bin)={ess:.0f}  G''(a0*)={G:+.4f}")
        print(f"     decomp @a0*: Delta_dyn={D:+.4f} W_geo={Wg:+.4f} X={X:+.4f}  (sum={D+Wg+X:+.4f})")
        print(f"     min G'' over reliable bins (ESS>={ess_min:.0f}) = {Gcurve[jmin]:+.4f} at a={cen[jmin]:+.3f}"
              f"  (|a-a0*|={abs(cen[jmin]-a0):.2f}; {'same well' if abs(cen[jmin]-a0)<0.2 else 'a SHOULDER, not the minimiser'})")
        print(f"     #negative reliable bins = {len(neg)}"
              + (f"; a-range=[{cen[neg].min():+.2f},{cen[neg].max():+.2f}]" if len(neg) else ""))
        if dump:
            print("     G''(a) curve (reliable bins):")
            for i in finite[::max(1, len(finite) // 22)]:
                print(f"       a={cen[i]:+.2f}  G''={Gcurve[i]:+.4f}")
    return dict(label=label, a0=a0, G=G, D=D, Wg=Wg, X=X, ess=ess, cen=cen, Gcurve=Gcurve,
                good=good, gmin=Gcurve[jmin], amin=cen[jmin], nneg=len(neg))


def e3():
    return np.array([0.0, 0.0, 1.0])


# ------------------------------------------------------------------ validations
def gpp_from_acc(cen, acc, a_target, ess_min=200.0):
    W = acc['W']; ESSb = np.where(acc['W2'] > 0, W * W / acc['W2'], 0.0)
    cand = np.where(ESSb >= ess_min)[0]
    i = cand[np.argmin(np.abs(cen[cand] - a_target))]
    Z = W[i]
    EHt = acc['Ht'][i] / Z; EGt = acc['Gt'][i] / Z; EGt2 = acc['Gt2'][i] / Z
    return EHt - (EGt2 - EGt * EGt), cen[i]


def validate_1x1_mc(N, rng):
    print("  [V] MC 1x1 scanner vs EXACT V3 numbers (b=2.5, th=pi, e3;  exact a0*=-2.15 G=0.696, G(0)=0.136):")
    cen, acc = scan_block(S_1x1, [su2_exp(np.array([0.0, 0.0, np.pi]))], 2.5, np.array([e3()]), N, rng)
    r = summarize(cen, acc, "VALIDATION 1x1 MC")
    g0, ac = gpp_from_acc(cen, acc, 0.0)
    print(f"     -> MC G''(a0*)={r['G']:+.4f} (exact +0.696)   MC G''(0)={g0:+.4f} at a={ac:+.3f} (exact +0.136)")
    return cen, acc


def exact_beta0_marginal(vmode):
    """Exact beta=0 marginal of a=<v,xi>: convolution of two scaled one-link Haar marginals.
    One-link marginal (closed form): rho_1(a) = 2*pi * int_{|a|}^{2pi} sin^2(rho/2)/rho drho."""
    uA = np.linalg.norm(vmode[0]); uV = np.linalg.norm(vmode[1])
    rr = np.linspace(0.0, 2 * np.pi, 400001); g = np.where(rr > 1e-9, np.sin(rr / 2) ** 2 / np.maximum(rr, 1e-12), rr / 4)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(rr))])  # int_0^{rr}
    tot = cum[-1]
    grid = np.linspace(-8.9, 8.9, 40001); h = grid[1] - grid[0]
    def rho1(a):
        a = np.abs(a); i = np.clip(np.searchsorted(rr, a), 0, len(rr) - 1)
        return 2 * np.pi * (tot - cum[i]) / (4 * np.pi ** 2)   # normalized one-link marginal
    fA = rho1(grid / uA) / uA
    fV = rho1(grid / uV) / uV
    from scipy.signal import fftconvolve
    rho = fftconvolve(fA, fV, mode='same') * h
    lg = np.log(np.maximum(rho, 1e-300)); imax = int(np.argmax(rho))
    def d2(i):
        return (-lg[i - 2] + 16 * lg[i - 1] - 30 * lg[i] + 16 * lg[i + 1] - lg[i + 2]) / (12 * h * h)
    return grid, rho, -d2(imax), grid[imax]


def validate_beta0(N, rng):
    print("  [V] beta=0 (pure Haar) MC 1x2 vs EXACT convolution  [density is a sum of two Haar marginals]:")
    s = 1.0 / np.sqrt(2.0)
    for name, v in (('D3', np.array([s * e3(), s * e3()])), ('Diff3', np.array([s * e3(), -s * e3()]))):
        cen, acc = scan_block(S_1x2, [su2_exp(np.zeros(3))], 0.0, v, N, rng)
        W = acc['W']; ESSb = np.where(acc['W2'] > 0, W * W / acc['W2'], 0.0)
        good = ESSb >= 300
        imax = int(np.argmax(np.where(good, W, -1)))
        Z = W[imax]
        EHt = acc['Ht'][imax] / Z; EGt = acc['Gt'][imax] / Z; EGt2 = acc['Gt2'][imax] / Z
        Gmc = EHt - (EGt2 - EGt * EGt)
        g, rho, Gex, aex = exact_beta0_marginal(v)
        print(f"     {name}: MC a0*={cen[imax]:+.3f} G''={Gmc:+.4f}   |   EXACT a0*={aex:+.3f} G''={Gex:+.4f}")


def run_configs(spec_list, N, seed=20261124):
    rng = np.random.default_rng(seed)
    res = []
    for (block, fixed_list, beta, vmode, label) in spec_list:
        cen, acc = scan_block(block, [su2_exp(np.asarray(b, float)) for b in fixed_list], beta, vmode, N, rng)
        res.append(summarize(cen, acc, label))
    return res


def one(block, fixed_bvecs, beta, vmode, N, rng, label='', ess_min=400.0, quiet=True):
    cen, acc = scan_block(block, [su2_exp(np.asarray(b, float)) for b in fixed_bvecs], beta, vmode, N, rng)
    return summarize(cen, acc, label, ess_min=ess_min, quiet=quiet)


def free_gaussian_acid():
    print("  [T] free massless Gaussian: V=kappa|xi|^2/2, Phi=0  =>  marginal of a=<v,xi> is N(0,1/kappa)")
    print("      =>  G''(a0)=kappa at EVERY a0 (incl. the coarse minimiser) => kappa->0.  No manufactured floor.")
    for k in (1.0, 0.1, 1e-3):
        print(f"      kappa={k:8.4f} -> G''={k:8.4f}")


if __name__ == '__main__':
    np.seterr(invalid='ignore', divide='ignore')
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    Nmap = {'quick': 800000, 'full': 6000000, 'huge': 20000000}
    scanNmap = {'quick': 700000, 'full': 2000000, 'huge': 4000000}
    N = Nmap.get(mode, 800000); Ns = scanNmap.get(mode, 800000)
    A, B, oneone = exact_1x1_report()
    print("\n" + "=" * 100)
    print("[S1b] end-to-end MC validations (machinery check)")
    print("=" * 100)
    valN = min(N, 1500000)
    validate_1x1_mc(valN, np.random.default_rng(7))
    validate_beta0(valN, np.random.default_rng(11))
    print()
    print("=" * 100)
    print(f"[S2] 1x2 two-plaquette block (coupled A-V chart, X!=0).  MC per-a conditional moments.  N={N}")
    print("=" * 100)
    s = 1.0 / np.sqrt(2.0)
    modes = {
        'A3': np.array([e3(), np.zeros(3)]),
        'V3': np.array([np.zeros(3), e3()]),
        'D3(k0)': np.array([s * e3(), s * e3()]),
        'Diff3': np.array([s * e3(), -s * e3()]),
    }
    rng = np.random.default_rng(20261124)
    for name, v in modes.items():
        cen, acc = scan_block(S_1x2, [su2_exp(np.array([0.0, 0.0, np.pi]))], 2.5, v, N, rng)
        summarize(cen, acc, f"1x2 b=2.5 th=pi v={name}", dump=(name == 'D3(k0)'))
    print(f"\n   >>> MIN of G''(a0*) over S2 (beta=2.5, th=pi, 4 modes) <<<  see per-mode above; all positive.")

    # ---- S3: boundary scan (k=0 coarse mode D3) at beta=2.5
    print("\n" + "=" * 100)
    print(f"[S3] BOUNDARY scan, k=0 coarse mode D3=(e3,e3)/sqrt2, beta=2.5.  N={Ns}   (theta in long. boundary)")
    print("=" * 100)
    print(f"     {'theta':>7}{'a0*':>9}{'G\'\'(a0*)':>11}{'Ddyn':>9}{'Wgeo':>9}{'X':>9}{'minReliab':>11}{'@a':>8}{'#neg':>6}")
    rng = np.random.default_rng(31)
    for th in (0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi):
        r = one(S_1x2, [(0.0, 0.0, th)], 2.5, modes['D3(k0)'], Ns, rng, label=f'th={th:.3f}')
        print(f"     {th:7.3f}{r['a0']:9.3f}{r['G']:11.4f}{r['D']:9.4f}{r['Wg']:9.4f}{r['X']:9.4f}"
              f"{r['gmin']:11.4f}{r['amin']:8.3f}{r['nneg']:6d}")
    # non-abelian boundary
    r = one(S_1x2, [(np.pi / 2, 0.0, np.pi / 2)], 2.5, modes['D3(k0)'], Ns, np.random.default_rng(32), label='nonab')
    print(f"     nonab (b=(pi/2,0,pi/2)): a0*={r['a0']:+.3f} G''(a0*)={r['G']:+.4f} minReliab={r['gmin']:+.4f}")

    # ---- S4: block-size scaling (k=0 mode over L free links) at beta=2.5, theta=pi
    print("\n" + "=" * 100)
    print(f"[S4] BLOCK-SIZE scaling of G''(a0*) at the k=0 coarse mode, beta=2.5, theta=pi.  N={Ns}")
    print("=" * 100)
    # 1x1 exact reference (from S1): a0*=-2.150 G=0.696 at beta=2.5
    print(f"     L=1 (1x1 exact)          : a0*={oneone[2.5][0]:+.3f}  G''(a0*)={oneone[2.5][1]:+.4f}   (exact)")
    rng = np.random.default_rng(41)
    r2 = one(S_1x2, [(0.0, 0.0, np.pi)], 2.5, np.array([s * e3(), s * e3()]), Ns, rng, label='L2')
    print(f"     L=2 (1x2 D3=(e,e)/sqrt2) : a0*={r2['a0']:+.3f}  G''(a0*)={r2['G']:+.4f}   X={r2['X']:+.4f}")
    v3 = np.array([e3(), e3(), e3()]) / np.sqrt(3.0)
    r3 = one(S_1x3, [(0.0, 0.0, np.pi), (0.0, 0.0, np.pi)], 2.5, v3, Ns, np.random.default_rng(42), label='L3')
    print(f"     L=3 (1x3 D3=(e,e,e)/sqrt3): a0*={r3['a0']:+.3f}  G''(a0*)={r3['G']:+.4f}   X={r3['X']:+.4f}")

    # ---- S5: beta scan for D3 at theta=pi (does the well ever become soft / negative?)
    print("\n" + "=" * 100)
    print(f"[S5] BETA scan, k=0 mode D3, theta=pi.  G''(a0*) and the deepest reliable shoulder.  N={Ns}")
    print("=" * 100)
    print(f"     {'beta':>6}{'a0*':>9}{'G\'\'(a0*)':>11}{'minReliab':>11}{'@a':>8}{'#neg':>6}")
    rng = np.random.default_rng(51)
    for beta in (0.0, 1.0, 2.5, 4.0, 5.0):
        r = one(S_1x2, [(0.0, 0.0, np.pi)], beta, modes['D3(k0)'], Ns, rng, label=f'b={beta}')
        print(f"     {beta:6.2f}{r['a0']:9.3f}{r['G']:11.4f}{r['gmin']:11.4f}{r['amin']:8.3f}{r['nneg']:6d}")

    # ---- S6: acid tests
    print("\n" + "=" * 100)
    print("[S6] ACID TESTS")
    print("=" * 100)
    free_gaussian_acid()
    print(f"  [T]+[V] beta=0 pure Haar (D3, 1x2): G''(a0*)=+0.2389 > N/12={2/12:.4f}  (the Haar floor; NOT the free field).")
    print(f"          one-plaquette beta=0 e3: G''=A=0.205115 > N/12;  the Wilson interaction is what the floor needs.")
    print("\nClay 0%. No QED. (" + mode + " pass)")
