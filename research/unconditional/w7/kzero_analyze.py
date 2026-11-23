"""
KZERO analysis: per-mode and pooled G''(0) with the (E[H]-Var)+W decomposition,
using block bootstrap over configs (robust to autocorrelation) for all SEs.

G''(0) KDE = curvature of -log P(a)  (gold standard, independent of identity).
decomposition at a=0 via Gaussian-kernel conditioning:
  bracket = E[H_vv|0]-Var(G_v|0),  W = E[Phi''|0]-Var(Phi'|0)-2Cov(G,Phi'|0)
  tot     = bracket + W   (must close to the KDE Gpp)
"""
import sys
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import rerun_core as RC  # noqa: E402


def _kde_point(a, bw_method=0.2, rfit=0.8):
    """KDE curvature of -log P(a) at 0; no internal bootstrap (outer bootstrap
    supplies the SE).  Same estimator as RC.kappa_kde's point value."""
    from scipy.stats import gaussian_kde
    a = np.asarray(a, float)
    sd = a.std()
    kde = gaussian_kde(a, bw_method=bw_method)
    xs = np.linspace(-rfit * sd, rfit * sd, 48)
    return float(2.0 * np.polyfit(xs, -np.log(kde(xs)), 2)[0])


def _kde(a, rng):
    return _kde_point(a)


def _decomp(a, H, G, P1, P2, bw):
    k = RC.kernel_cond(a, {'H': H, 'G': G, 'P1': P1, 'P2': P2}, bw=bw)
    EH, VG = k['E_H'], k['V_own']['G']
    EP2, VP1, cov = k['E_P2'], k['V_own']['P1'], k['C_G_P1']
    return dict(EH=EH, VG=VG, EP2=EP2, VP1=VP1, cov=cov,
                bracket=EH - VG, W=EP2 - VP1 - 2.0 * cov,
                tot=(EH + EP2) - (VG + VP1 + 2.0 * cov))


def block_bootstrap(N, nboot, block, rng):
    """yield arrays of resampled indices, contiguous blocks (circular)."""
    nb = int(np.ceil(N / block))
    for _ in range(nboot):
        starts = rng.integers(0, N, nb)
        idx = np.concatenate([(np.arange(s, s + block) % N) for s in starts])[:N]
        yield idx


def ac_time(x):
    """integrated autocorrelation time (units of stored configs)."""
    x = x - x.mean()
    v = np.dot(x, x) / x.size
    if v == 0:
        return 1.0
    s = 1.0
    for lag in range(1, min(200, x.size // 4)):
        c = np.dot(x[:-lag], x[lag:]) / (x.size - lag) / v
        if c <= 0:
            break
        s += 2.0 * c
    return max(s, 1.0)


def analyse(path, bw=0.30, nboot=300, block=20, rngseed=7000, verbose=True):
    d = np.load(path)
    A, H, G, P1, P2 = d['a'], d['H'], d['G'], d['P1'], d['P2']
    N, Nm = A.shape
    names = [str(s) for s in d['names']] if 'names' in d else [f'm{j}' for j in range(Nm)]
    L = int(d['L']) if 'L' in d else -1
    if 'kvals' in d:
        kvals = np.asarray(d['kvals'], float)
    else:
        kphys = float(d['kphys']) if 'kphys' in d else (2 * np.pi / L if L > 0 else np.nan)
        kvals = np.full(Nm, kphys)
    # group mode indices by momentum
    uk = np.unique(np.round(kvals, 6))
    rng = np.random.default_rng(rngseed)
    tau = np.array([ac_time(A[:, j]) for j in range(Nm)])
    rows = []
    for j in range(Nm):
        a, h, g, p1, p2 = A[:, j], H[:, j], G[:, j], P1[:, j], P2[:, j]
        kap = _kde(a, rng)
        dec = _decomp(a, h, g, p1, p2, bw)
        # block bootstrap
        bk = {k: [] for k in ('kap', 'bracket', 'W', 'tot', 'EH', 'VG', 'EP2', 'VP1', 'cov')}
        for idx in block_bootstrap(N, nboot, block, rng):
            aa, hh, gg, q1, q2 = a[idx], h[idx], g[idx], p1[idx], p2[idx]
            try:
                bk['kap'].append(_kde(aa, rng))
            except Exception:
                continue
            dc = _decomp(aa, hh, gg, q1, q2, bw)
            for k in ('bracket', 'W', 'tot', 'EH', 'VG', 'EP2', 'VP1', 'cov'):
                bk[k].append(dc[k])
        se = {k: float(np.std(v)) for k, v in bk.items()}
        rows.append(dict(name=names[j], kap=kap, se=se['kap'], inv=1.0 / np.var(a),
                         sd=float(a.std()), **dec,
                         s_bracket=se['bracket'], s_W=se['W'], s_tot=se['tot'],
                         s_EH=se['EH'], s_VG=se['VG'], s_EP2=se['EP2'],
                         s_VP1=se['VP1'], s_cov=se['cov'], tau=tau[j]))
    # pooled estimate (equal weight over modes)
    def pool(key):
        v = np.array([r[key] for r in rows])
        return float(v.mean()), v
    kap_m, kaps = pool('kap')
    # pooled SE = combine mean bootstrap SE with mode-to-mode scatter
    se_boot = float(np.sqrt(np.mean([r['se'] ** 2 for r in rows]) / Nm))
    se_mode = float(kaps.std(ddof=1) / np.sqrt(Nm))
    se_pool = float(np.sqrt(se_boot ** 2 + se_mode ** 2)) if Nm > 1 else se_boot
    pooled = dict(kap=kap_m, se=se_pool, se_boot=se_boot, se_mode=se_mode)
    for key in ('bracket', 'W', 'tot', 'EH', 'VG', 'EP2', 'VP1', 'cov', 'inv'):
        m, _ = pool(key)
        pooled[key] = m
    if verbose:
        print(f"# {path}  L={L} beta={d['beta'] if 'beta' in d else '?'} "
              f"N={N} modes={Nm}")
        print(f"  {'mode':18s} {'k':>6s} {'Gpp':>8s} {'+-SE':>7s} {'1/Var':>6s} "
              f"{'brkt':>7s} {'W':>7s} {'tot':>7s} {'E[H]':>6s} {'VarG':>6s} "
              f"{'EPhi2':>6s} {'tau':>5s}")
        for j, r in enumerate(rows):
            print(f"  {r['name']:18s} {kvals[j]:6.3f} {r['kap']:+8.4f} {r['se']:7.4f} {r['inv']:6.3f} "
                  f"{r['bracket']:+7.3f} {r['W']:+7.3f} {r['tot']:+7.3f} "
                  f"{r['EH']:6.3f} {r['VG']:6.3f} {r['EP2']:6.3f} {r['tau']:5.2f}")
        print(f"  POOLED  Gpp={pooled['kap']:+.4f} +- {pooled['se']:.4f} "
              f"(boot {pooled['se_boot']:.4f} + mode {pooled['se_mode']:.4f})  "
              f"1/Var={pooled['inv']:.3f}")
        print(f"          bracket={pooled['bracket']:+.3f}  W={pooled['W']:+.3f}  "
              f"tot={pooled['tot']:+.3f}  E[H]={pooled['EH']:.3f} "
              f"VarG={pooled['VG']:.3f} EPhi2={pooled['EP2']:.3f} "
              f"(closure tot-Gpp={pooled['tot']-pooled['kap']:+.3f})")
        # per-momentum pooled (for k-scans on one volume)
        if len(uk) > 1:
            print(f"  --- per-k ---")
            for kv in uk:
                js = [j for j in range(Nm) if abs(kvals[j] - kv) < 1e-6]
                kk = np.array([rows[j]['kap'] for j in js])
                se_b = np.sqrt(np.mean([rows[j]['se'] ** 2 for j in js]) / len(js))
                se_m = kk.std(ddof=1) / np.sqrt(len(js)) if len(js) > 1 else 0.0
                bm = np.mean([rows[j]['bracket'] for j in js])
                wm = np.mean([rows[j]['W'] for j in js])
                tm = np.mean([rows[j]['tot'] for j in js])
                eh = np.mean([rows[j]['EH'] for j in js])
                vg = np.mean([rows[j]['VG'] for j in js])
                ep = np.mean([rows[j]['EP2'] for j in js])
                print(f"  k={kv:6.3f} ({len(js)} modes) Gpp={kk.mean():+.4f} +- "
                      f"{np.sqrt(se_b**2 + se_m**2):.4f}  1/Var={np.mean([rows[j]['inv'] for j in js]):.3f}  "
                      f"brkt={bm:+.3f} W={wm:+.3f} tot={tm:+.3f}  E[H]={eh:.3f} VarG={vg:.3f} EPhi2={ep:.3f}")
    return rows, pooled


if __name__ == '__main__':
    import argparse
    import os
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='+')
    ap.add_argument('--block', type=int, default=20)
    ap.add_argument('--nboot', type=int, default=300)
    ap.add_argument('--bw', type=float, default=0.30)
    a = ap.parse_args()
    for fn in a.files:
        analyse(fn, bw=a.bw, nboot=a.nboot, block=a.block)
        print()
