"""
H1-NUMERIC-RIGHT: measure G'' on the CORRECT object, not the gauge-dependent gluon.

WHY THIS SCRIPT EXISTS (read K0-OVERLAP-AUDIT.md + NUMERICS-AUDIT.md first)
---------------------------------------------------------------------------
The campaign's +0.3 curvature (COMPUTE-GPP-K0-CONTINUUM) was measured on
    a = < v , log U >     (v = one transverse plane wave, one polarisation, one colour)
This is a gauge-DEPENDENT, colour-ADJOINT, 1^- VECTOR (a single gluon amplitude),
NOT the gauge-invariant 0^{++} scalar.  Its G'' is carried by the Haar/Jacobian
term W_geo, which the FREE theory satisfies too and which is LARGER for the free
field at the same k.  So the +0.3 never engages the acid test (free -> 0 only at k->0).

This script measures TWO things, both honestly labelled:

  (A) GLUON  (the cone/vector object; H1-relevant per K0-OVERLAP-AUDIT (c)(2),
      because H1 is a measure-level LSI that INCLUDES the cone, but NOT the 0^{++}).
      Across a beta grid {1.0,1.5,2.0,2.2,2.4} SU(2) L=8 we report, per beta:
        G''(KDE), bracket=E[H]-Var(G), W_geo=E[Phi'']-Var(Phi'), X=-2Cov(G,Phi'),
        S=1/Var(a), D/S = 1 - G''*Var  (the non-Gaussianity defect; W7 anchor ~2%),
        free value 4 sin^2(k/2) at the SAME k (the acid test side by side), and
        one-sided lower bounds.  This is the table the request literally asked for,
        on the only object where the exact identity holds -- clearly caveated.

  (B) 0^{++} DENSITY  (the request's suggested "correct" object).
      q(x) = 1 - (1/6) sum_{mu<nu} (1/2) Re tr U_p(x)  ~  tr F^2, gauge-invariant.
      We record the k=0 mode (qbar = volume mean) and the softest nonzero Fourier
      amplitude q~(k1) (real cos quadrature).  KEY POINT, tested numerically:
        Var(qbar) = chi / V,  Var(q~(k)) = C~(k) / (2V)     (chi, C~ = susceptibilities)
      so the gauge-invariant density's amplitude variance VANISHES as V -> 0^{-1}? no,
      as V -> infinity (like 1/V), hence its "G'' = 1/Var" DIVERGES with volume.
      That is a SUSCEPTIBILITY, the OPPOSITE of an H1 amplitude (which needs Var not->0).
      We demonstrate this with a volume scan L = 4,6,8 at fixed beta.

      We also extract the 0^{++} glueball mass m0 from the connected correlator
      <q(x)q(0)>_c -- the genuine gauge-invariant 0^{++} gap -- to show the 0^{++}
      IS gapped, but that this is a 2-point-function (decay-rate) statement, not a
      G''-curvature sign.  Exactly the K0-OVERLAP-AUDIT (c)(3) caveat.

MEMORY DISCIPLINE (charter sec.4): ONE process, ONE lattice; configs NEVER stored
(streaming in-place heatbath, single U holder ~0.5 MB); only N x Nmode scalars and
an L-vector correlator are accumulated.  Jobs run SEQUENTIALLY (never parallel).
Peak RSS is reported per job.  L <= 8, SU(2).
"""
import sys
import os
import time
import json
import ctypes
import ctypes.wintypes as wt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
W7 = 'D:/Dev/math/research/unconditional/w7'
CPB = 'D:/Dev/math/research/lattice-cpb'
sys.path.insert(0, CPB)
sys.path.insert(0, W7)
sys.path.insert(0, HERE)

import w4_blocked as W          # noqa: E402
import gpp_core as C            # noqa: E402
import rerun_core as RC         # noqa: E402
import su3 as S3               # noqa: E402  (su2_heatbath lives here)
import k0_continuum_analyze as KA   # noqa: E402  (stats, delete_block_jack, pooled_jack)
import corr_su2 as CS           # noqa: E402  (plaq_field, corr, fit_m)

adj = W.adj

# 6 transverse plane-wave modes (one polarisation, one colour each), SU(2)
SPEC2 = [(1, 2, 2), (1, 3, 2), (0, 2, 2), (2, 3, 2), (3, 1, 2), (1, 2, 0)]


# ----------------------------------------------------------------------
# peak RSS (Windows), MB
# ----------------------------------------------------------------------
class _PMC(ctypes.Structure):
    _fields_ = [('cb', wt.DWORD), ('PageFaultCount', wt.DWORD),
                ('PeakWorkingSetSize', ctypes.c_size_t),
                ('WorkingSetSize', ctypes.c_size_t),
                ('QuotaPeakPagedPoolUsage', ctypes.c_size_t),
                ('QuotaPagedPoolUsage', ctypes.c_size_t),
                ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t),
                ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                ('PagefileUsage', ctypes.c_size_t),
                ('PeakPagefileUsage', ctypes.c_size_t)]


_PSAPI = ctypes.WinDLL('psapi', use_last_error=True)
_PSAPI.GetProcessMemoryInfo.argtypes = [wt.HANDLE, ctypes.POINTER(_PMC), wt.DWORD]
_PSAPI.GetProcessMemoryInfo.restype = wt.BOOL
_K32 = ctypes.WinDLL('kernel32', use_last_error=True)
_K32.GetCurrentProcess.restype = wt.HANDLE


def _pmc():
    pmc = _PMC(); pmc.cb = ctypes.sizeof(_PMC)
    _PSAPI.GetProcessMemoryInfo(_K32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb)
    return pmc


def peak_rss_mb():
    return _pmc().PeakWorkingSetSize / 1024.0 / 1024.0


def rss_mb():
    return _pmc().WorkingSetSize / 1024.0 / 1024.0


# ----------------------------------------------------------------------
# streaming SU(2) sampler (in-place heatbath, single config holder)
# ----------------------------------------------------------------------
def make_sampler(L, beta, seed, cold=False):
    rng = np.random.default_rng(seed)
    if cold:
        U = np.zeros((4, L, L, L, L, 2, 2), dtype=np.complex128)
        for mu in range(4):
            U[mu, ..., 0, 0] = 1.0
            U[mu, ..., 1, 1] = 1.0
    else:
        U = W.random_su2(rng, (4, L, L, L, L))
    ns = L ** 4
    x0, x1, x2, x3 = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
    par = (x0 + x1 + x2 + x3) % 2
    idxp = [np.nonzero(par.reshape(-1) == c)[0] for c in (0, 1)]
    kappa = beta / 2.0

    def sweep():
        Uf = U.reshape(4, ns, 2, 2)
        for c in (0, 1):
            ids = idxp[c]
            for mu in range(4):
                S = W.staple_sum2(U, mu).reshape(ns, 2, 2)[ids]
                A = Uf[mu, ids]
                M = A @ S
                v = S3.su2_heatbath(kappa, M, 0, 1, rng)
                Uf[mu, ids] = v @ A

    return U, sweep


# ----------------------------------------------------------------------
# per-config measurement:  gluon decomposition  +  gauge-invariant 0^{++} density
# ----------------------------------------------------------------------
def build_gluon_modes(L, nmax):
    modes, names, kvals = [], [], []
    for n in range(1, nmax + 1):
        for (ax, pl, c) in SPEC2:
            modes.append(C.plane_mode(L, ax, pl, c, n=n, kind='cos'))
            names.append(f'ax{ax}_p{pl}_c{c}_n{n}')
            kvals.append(2.0 * np.pi * n / L)
    return modes, names, np.array(kvals)


def measure(U, beta, gmodes, L, nmax_dens=2):
    # ---- gluon (flat chart X=log U, a=<v,x>) ----
    x = C.su2_log_x(U)
    Xfrom = RC.X_from_x_su2
    S0 = W.action2(U, beta)
    a = []; H = []; G = []; P1 = []; P2 = []
    for v in gmodes:
        a.append(float(np.sum(v * x)))
        Sp = W.action2(RC.expi(Xfrom(x + 1e-3 * v)), beta)
        Sm = W.action2(RC.expi(Xfrom(x - 1e-3 * v)), beta)
        G.append((Sp - Sm) / 2e-3)
        H.append((Sp - 2.0 * S0 + Sm) / 1e-6)
        p1, p2 = RC.phi_derivs_analytic(x, v, Xfrom)
        P1.append(p1); P2.append(p2)

    # ---- gauge-invariant 0^{++} density s(x) = (1/6) sum 1/2 Re tr U_p ----
    s = CS.plaq_field(U)                       # (L,L,L,L), gauge invariant
    V = L ** 4
    qbar = float(s.mean())                     # k=0 mode of the density
    # softest nonzero Fourier amplitudes (cos quadrature), per axis, n=1..nmax_dens
    amps = {}                                  # key (n,ax) -> (Ac, As)
    xs = np.arange(L)
    for n in range(1, nmax_dens + 1):
        k = 2.0 * np.pi * n / L
        for ax in range(3):
            shp = [L if i == ax else 1 for i in range(4)]
            ph = (k * xs).reshape(shp)
            cc = np.broadcast_to(np.cos(ph), s.shape)
            ss = np.broadcast_to(np.sin(ph), s.shape)
            amps[(n, ax)] = (float(np.sum(s * cc)) / V, float(np.sum(s * ss)) / V)
    # connected correlator along axis 0
    cr = CS.corr(s, axis=0)                    # (L,)

    return dict(a=np.array(a), H=np.array(H), G=np.array(G),
                P1=np.array(P1), P2=np.array(P2),
                qbar=qbar, amps=amps, cr=cr)


# ----------------------------------------------------------------------
# density statistics with delete-block jackknife
# ----------------------------------------------------------------------
def dens_stats(cfg, L):
    """cfg: list of per-config dicts.  Returns full stats + jackknife SE."""
    N = len(cfg)
    qbar = np.array([c['qbar'] for c in cfg])
    V = L ** 4
    # pooled nonzero amplitudes at n=1 (3 axes) -> marginal sample
    A1 = np.concatenate([np.array([c['amps'][(1, ax)][0] for c in cfg])
                         for ax in range(3)])     # (3N,)
    cr = np.mean([c['cr'] for c in cfg], axis=0)

    def stats_of(qb, A):
        varA = float(np.var(A))
        Gpp, _ = RC.kappa_kde(A, rng=np.random.default_rng(0), nboot=60, bw_method=0.2)
        S = 1.0 / varA
        return dict(Var_qbar=float(np.var(qb)), chi=float(np.var(qb)) * V,
                    Var_A=varA, Gpp_dens=float(Gpp), S=S,
                    DS=float(1.0 - Gpp * varA))

    full = stats_of(qbar, A1)
    # jackknife over configs (delete 1/nblk of the configs; axes travel with config)
    nblk = int(min(20, max(4, N // 8)))
    edges = np.linspace(0, N, nblk + 1).astype(int)
    A_by_ax = [np.array([c['amps'][(1, ax)][0] for c in cfg]) for ax in range(3)]
    leave = {k: [] for k in full}
    for b in range(nblk):
        keep = np.ones(N, bool); keep[edges[b]:edges[b + 1]] = False
        Ab = np.concatenate([A_by_ax[ax][keep] for ax in range(3)])
        st = stats_of(qbar[keep], Ab)
        for k in full:
            leave[k].append(st[k])
    se = {}
    for k in full:
        arr = np.array(leave[k])
        se[k] = float(np.sqrt((nblk - 1) / nblk * np.sum((arr - arr.mean()) ** 2)))
    # glueball mass from the connected correlator
    try:
        m0 = CS.fit_m(cr, L)
    except Exception:
        m0 = float('nan')
    meff1 = float('nan')
    if min(cr[0], cr[1], cr[2]) > 0:
        meff1 = float(np.arccosh((cr[0] + cr[2]) / (2 * cr[1])))
    return dict(full=full, se=se, cr=cr, m0=m0, meff1=meff1, nblk=nblk)


# ----------------------------------------------------------------------
# one beta: stream N configs, collect gluon + density
# ----------------------------------------------------------------------
def run_beta(L, beta, N, gap, therm, seed, nmax_gluon=1, nmax_dens=2, want_density=True):
    gmodes, gnames, kvals = build_gluon_modes(L, nmax_gluon)
    Nm = len(gmodes)
    print(f"# H1-RIGHT {beta=} L={L} N={N} gap={gap} therm={therm} seed={seed} "
          f"gluon_modes={Nm}", flush=True)
    U, sweep = make_sampler(L, beta, seed, cold=False)
    t0 = time.time()
    for _ in range(therm):
        sweep()
    plaq = float(1.0 - W.action2(U, beta) / (beta * 6.0 * L ** 4))
    print(f"  therm done {time.time()-t0:.0f}s plaq={plaq:.4f} "
          f"RSS={rss_mb():.0f} peak={peak_rss_mb():.0f}MB", flush=True)

    acc = {k: np.empty((N, Nm)) for k in ('a', 'H', 'G', 'P1', 'P2')}
    cfg = [] if want_density else None
    t0 = time.time(); tlast = t0
    for i in range(N):
        for _ in range(gap):
            sweep()
        r = measure(U, beta, gmodes, L, nmax_dens=nmax_dens)
        for k in acc:
            acc[k][i] = r[k]
        if want_density:
            cfg.append(dict(qbar=r['qbar'], amps=r['amps'], cr=r['cr']))
        if (i + 1) % 50 == 0 or i == 0:
            dt = time.time() - tlast; tlast = time.time()
            eta = (N - i - 1) / (i + 1) * (time.time() - t0)
            print(f"    {i+1}/{N} ({time.time()-t0:.0f}s ETA {eta:.0f}s) "
                  f"RSS={rss_mb():.0f} peak={peak_rss_mb():.0f}MB", flush=True)

    out = dict(beta=beta, L=L, N=N, nmodes=Nm, plaq=plaq,
               peak_rss_mb=peak_rss_mb(), gluon=None)
    if Nm > 0:
        # gluon decomposition (reuse validated jackknife)
        full, se, ml, nblk = KA.pooled_jack(acc['a'], acc['H'], acc['G'],
                                            acc['P1'], acc['P2'], nblk=20)
        k0 = float(np.mean(kvals))
        free = 4.0 * np.sin(0.5 * k0) ** 2        # free transverse G'' at k=2pi/L
        inv = full['inv']                         # inv = S = 1/Var(a)
        DS = float(1.0 - full['Gpp'] / inv)       # D/S = 1 - G''*Var = 1 - G''/S
        lb_stat = float(full['Gpp'] - 1.645 * se['Gpp'])
        lb_rig = float(full['brX'] - 1.645 * se['brX'])
        out['k'] = k0
        out['gluon'] = dict(full=full, se=se, DS=DS, free=free,
                            ratio_Gpp_over_free=float(full['Gpp'] / free),
                            LB_stat=lb_stat, LB_rig=lb_rig)
        print(f"  GLUON beta={beta}: Gpp={full['Gpp']:+.4f}+-{se['Gpp']:.4f} "
              f"({full['Gpp']/se['Gpp']:+.1f}s) bracket={full['bracket']:+.3f} "
              f"Wgeo={full['Wgeo']:+.3f} X={full['X']:+.3f} S=1/Var={inv:.3f} "
              f"D/S={DS:.3f} free={free:.3f} G/free={full['Gpp']/free:+.2f} "
              f"LB_rig={lb_rig:+.3f}", flush=True)

    if want_density:
        d = dens_stats(cfg, L)
        out['density'] = d
        f = d['full']; s = d['se']
        print(f"  DENSITY beta={beta}: Var(qbar)={f['Var_qbar']:.3e} "
              f"chi=V*Var={f['chi']:.4f} | Var_A={f['Var_A']:.3e} "
              f"Gpp_dens={f['Gpp_dens']:+.3e} S=1/Var={f['S']:.3e} "
              f"D/S={f['DS']:+.3f} | m0={d['m0']:.3f} meff1={d['meff1']:.3f}",
              flush=True)
    print(f"  DONE beta={beta} {time.time()-t0:.0f}s peak RSS={peak_rss_mb():.0f}MB",
          flush=True)
    return out


# ----------------------------------------------------------------------
def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--betas', type=float, nargs='+',
                    default=[1.0, 1.5, 2.0, 2.2, 2.4])
    ap.add_argument('--L', type=int, default=8)
    ap.add_argument('--N', type=int, default=200)
    ap.add_argument('--gap', type=int, default=4)
    ap.add_argument('--therm', type=int, default=250)
    ap.add_argument('--seed', type=int, default=4242)
    ap.add_argument('--nmax_gluon', type=int, default=1)
    ap.add_argument('--vol_betas', type=float, nargs='*', default=[2.4, 1.5],
                    help='betas at which to run L=4,6 density-only for the 1/V test '
                         '(L=8 main run supplies the third point)')
    ap.add_argument('--Nvol', type=int, default=120)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()

    results = dict(meta=dict(L=a.L, N=a.N, gap=a.gap, therm=a.therm,
                             seed=a.seed, nmax_gluon=a.nmax_gluon),
                   per_beta=[], volume_scan=[])
    t_all = time.time()
    for beta in a.betas:
        results['per_beta'].append(run_beta(a.L, beta, a.N, a.gap, a.therm,
                                            a.seed + int(1000 * beta),
                                            nmax_gluon=a.nmax_gluon))
    # volume scan (L=4,6 density-only; L=8 from per_beta) for the Var~1/V test
    for beta in a.vol_betas:
        for LL in (4, 6):
            r = run_beta(LL, beta, a.Nvol, a.gap, a.therm,
                         a.seed + int(1000 * beta) + LL,
                         nmax_gluon=0, want_density=True)
            results['volume_scan'].append(dict(beta=beta, L=LL,
                                               density=r['density']))

    def _js(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.integer):
            return int(o)
        return float(o)

    out = a.out or os.path.join(HERE, f'h1_right_L{a.L}_N{a.N}.json')
    with open(out, 'w') as f:
        json.dump(results, f, indent=1, default=_js)
    print(f"\n# wrote {out}  total {time.time()-t_all:.0f}s "
          f"peak RSS={peak_rss_mb():.0f}MB", flush=True)


if __name__ == '__main__':
    main()
