"""
W7 core: EXACT marginal curvature G''(0) of the lowest transverse gluon mode.

Object (from MARGINAL.md / W6):
    G''(a) = E[ v^T H v | a ] - Var( v^T grad S | a )
with a = v^T x  the amplitude of ONE collective mode v (the k->0 / lowest
transverse direction), and the conditional average over everything else.

Exactness route (no Gaussian proxy):  define the 1-D marginal distribution
    P(a) = prob( v^T x = a )   from the equilibrium ensemble,
    V_eff(a) = -log P(a) + const.
Then EXACTLY (constrained free energy):
    G''(a) = V_eff''(a) = E[H_vv | a] - Var(G_v | a),
with H_vv = v^T H v (curvature of S along v) and G_v = v^T grad S (force).
So the equilibrium histogram of a gives G''(0) exactly (up to statistics),
with NO conditional-MCMC and NO linear-regression proxy.

We cross-check three ways per mode:
  (i)  kappa_hist  = curvature of -log P(a) near a=0 (gold standard);
  (ii) kappa_invvar= 1/Var(a)  (equals kappa_hist for a Gaussian mode);
  (iii) kappa_Ev    = E[H_vv] - Var(G_v)  (the literal estimator, unconditioned);
plus binned conditional <H_vv|bin> - Var(G_v|bin).

Conventions: SU(2) Wilson action S=beta*sum_p(1-1/2 Re tr U_p);
t^a=sigma^a/2, left variation U->exp(i xi_a t^a)U, x_a = 2 Im tr(t^a log U).
Modes are UNIT-normalised in the (4,L^4,3) real algebra space.
"""
import sys
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import w4_blocked as W  # noqa: E402

T2 = W.T2  # (3,2,2)


# ----------------------------------------------------------------------
# SU(2) principal-log algebra components  x_a :  U = exp(i x_a t^a)
# ----------------------------------------------------------------------
def su2_log_x(U):
    tr = np.real(np.trace(U, axis1=-2, axis2=-1))          # 2 cos theta
    c = np.clip(tr / 2.0, -1.0, 1.0)
    th = np.arccos(c)                                      # in [0, pi]
    s = np.sin(th)
    UmU = U - W.adj(U)
    trt = np.einsum('aij,...ji->...a', T2, UmU)            # complex
    fac = np.where(s > 1e-7, th / np.maximum(s, 1e-12), 1.0)
    x = fac[..., None] * trt.imag
    # singular points U ~ -I (s->0, theta->pi): set to 0 (measure ~0); flag count
    return x


# ----------------------------------------------------------------------
# exp(i * comp . t) @ U   (left multiplication along the ray)
# ----------------------------------------------------------------------
def apply_exp(U, comp):
    """U -> exp(i comp.t) U, comp: (4,L,L,L,L,3)."""
    A = np.einsum('...a,aij->...ij', comp, T2)             # Hermitian
    w, V = np.linalg.eigh(A)
    E = np.einsum('...ik,...k,...jk->...ij', V, np.exp(1j * w), np.conj(V))
    return E @ U


def ray_action(U, comp, beta, s):
    return W.action2(apply_exp(U, s * comp), beta)


# ----------------------------------------------------------------------
# plane-wave transverse mode  v  (unit-normalised, one pol, one colour)
# ----------------------------------------------------------------------
def plane_mode(L, axis, pol, color, n=1, kind='cos'):
    """v[mu,x,a]; mu in 0..3 (spacetime), pol != axis (transverse), one colour."""
    xs = np.arange(L)
    xa = xs.reshape([L if i == axis else 1 for i in range(4)])
    ph = (np.cos if kind == 'cos' else np.sin)(2.0 * np.pi * n * xa / L)
    ph = np.broadcast_to(ph, (L, L, L, L)).astype(np.float64)
    v = np.zeros((4, L, L, L, L, 3))
    v[pol, :, :, :, :, color] = ph
    v /= np.linalg.norm(v)
    return v


# ----------------------------------------------------------------------
# per-config observables for a set of modes
# ----------------------------------------------------------------------
def mode_obs(U, beta, modes, h=1e-3):
    """returns dict arrays: a[k], Hvv[k], Gv[k] for each mode."""
    x = su2_log_x(U)                       # (4,L,L,L,L,3)
    S0 = W.action2(U, beta)
    a, Hvv, Gv = [], [], []
    for v in modes:
        a.append(float(np.sum(v * x)))
        # force & curvature along the ray U(s)=exp(i s v)U
        Sp = ray_action(U, v, beta, h)
        Sm = ray_action(U, v, beta, -h)
        Gv.append((Sp - Sm) / (2.0 * h))
        Hvv.append((Sp - 2.0 * S0 + Sm) / (h * h))
    return np.array(a), np.array(Hvv), np.array(Gv)


# ----------------------------------------------------------------------
# V_eff''(0) from a scalar sample by even-polynomial fit of -log P(a)
# ----------------------------------------------------------------------
def veff_curvature(a, nboot=400, rng=None, nave=2.5, nbins=40):
    """Fit -log p(a) = c + 0.5*kappa*a^2 + (lam/4) a^4 on |a|<nave*sigma.
    Returns kappa, its bootstrap SE, and lam.  Also invvar=1/Var(a)."""
    a = np.asarray(a, dtype=float)
    sd = a.std()
    if rng is None:
        rng = np.random.default_rng(0)

    def fit(sample):
        m = np.abs(sample) < nave * sd
        if m.sum() < nbins:
            m = np.abs(sample) < 4 * sd
        lo, hi = sample[m].min(), sample[m].max()
        edges = np.linspace(lo, hi, nbins + 1)
        cnt, _ = np.histogram(sample[m], bins=edges, density=False)
        cen = 0.5 * (edges[1:] + edges[:-1])
        good = cnt > 0
        cc = cen[good]
        yy = -np.log(cnt[good] / (sample[m].size * (edges[1] - edges[0])))
        # even polynomial: yy = c0 + c1*a^2 + c2*a^4
        A = np.vstack([np.ones_like(cc), cc ** 2, cc ** 4]).T
        coef, *_ = np.linalg.lstsq(A, yy, rcond=None)
        return 2.0 * coef[1], 4.0 * coef[2]   # kappa, lam

    kap, lam = fit(a)
    boots = np.empty(nboot)
    n = a.size
    for b in range(nboot):
        idx = rng.integers(0, n, n)
        boots[b] = fit(a[idx])[0]
    return kap, boots.std(), lam, float(np.mean(boots)), float(1.0 / a.var())


def ev_estimator(a, Hvv, Gv, nboot=400, rng=None, nbins=12, nave=1.5):
    """E[Hvv]-Var(Gv) unconditioned + binned conditional near a=0."""
    a = np.asarray(a, float); Hvv = np.asarray(Hvv, float); Gv = np.asarray(Gv, float)
    if rng is None:
        rng = np.random.default_rng(1)
    Ev = Hvv.mean() - Gv.var()
    # bootstrap of Ev by resampling configs
    n = a.size
    bs = np.empty(nboot)
    for b in range(nboot):
        idx = rng.integers(0, n, n)
        bs[b] = Hvv[idx].mean() - Gv[idx].var()
    # conditional near a=0: bin, use bin containing 0
    sd = a.std()
    m0 = np.abs(a) < nave * sd / np.sqrt(nbins)
    if m0.sum() < 8:
        m0 = np.abs(a) < np.quantile(np.abs(a), 0.15)
    Ev0 = Hvv[m0].mean() - Gv[m0].var()
    return Ev, bs.std(), Ev0, float(Hvv.mean()), float(Gv.var()), int(m0.sum())


def sample_configs(L, beta, ntherm, ncfg, gapevery, seed, cold=True):
    rng = np.random.default_rng(seed)
    return W.sample2(L, beta, rng, ntherm, ncfg, gapevery)


if __name__ == '__main__':
    import time
    L = 4
    beta = 2.5
    ncfg = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    t0 = time.time()
    cfgs = sample_configs(L, beta, 400, ncfg, 3, seed=101)
    print(f"sampled {ncfg} cfgs in {time.time()-t0:.1f}s; plaq="
          f"{np.mean([W.action2(U,beta) for U in cfgs]):.3e}", flush=True)
    # a few transverse plane modes at the lowest momentum k=2pi/L
    modes = {
        'k1_pol2_c3_cos': plane_mode(L, axis=1, pol=2, color=2, n=1, kind='cos'),
        'k1_pol3_c3_cos': plane_mode(L, axis=1, pol=3, color=2, n=1, kind='cos'),
        'k1_pol0_c3_cos': plane_mode(L, axis=1, pol=0, color=2, n=1, kind='cos'),
        'k1_pol2_c1_cos': plane_mode(L, axis=1, pol=2, color=0, n=1, kind='cos'),
        'k1_pol2_c3_sin': plane_mode(L, axis=1, pol=2, color=2, n=1, kind='sin'),
    }
    names = list(modes)
    ml = [modes[n] for n in names]
    data = {n: ([], [], []) for n in names}
    for i, U in enumerate(cfgs):
        a, H, G = mode_obs(U, beta, ml, h=1e-3)
        for n, av, hv, gv in zip(names, a, H, G):
            data[n][0].append(av); data[n][1].append(hv); data[n][2].append(gv)
    rng = np.random.default_rng(7)
    print(f"\n{'mode':18s} {'kappa_hist':>11s} {'+-SE':>8s} {'1/Var':>8s} "
          f"{'E[H]':>8s} {'Var(G)':>8s} {'E[H]-Var':>9s}")
    for n in names:
        a, H, G = map(np.array, data[n])
        kap, se, lam, mb, inv = veff_curvature(a, rng=rng)
        Ev, Ese, Ev0, EH, VG, nb = ev_estimator(a, H, G, rng=rng)
        print(f"{n:18s} {kap:+11.4f} {se:8.4f} {inv:8.4f} "
              f"{EH:8.4f} {VG:8.4f} {Ev:+9.4f}   lam={lam:+.3f} nbin0={nb}")
