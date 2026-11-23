"""Pool + robustness analysis for the SU(3) k=2pi/L curvature computation.

Loads any su3_k2pil_L{L}_N*_s*.npz (fresh) plus the pre-existing
rerun_su3_L{L}.npz, pools the samples per (L, mode), and reports:
  - G'' via KDE curvature of -log P(a), with bootstrap SE (gold standard)
  - G'' via 1/Var(a)
  - the exact-identity decomposition at a=0: bracket=E[H]-Var(G), W, bracket+W
  - bandwidth sensitivity of the decomposition (bw=0.20/0.30/0.40)
  - lag-1 autocorrelation of the collective amplitude a (decorrelation check)
Pools fresh-only and fresh+existing separately.
"""
import sys
import glob
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import rerun_run as R          # noqa: E402
import rerun_core as RC        # noqa: E402

SPEC = [(1, 2, 3), (1, 3, 3), (0, 2, 3), (2, 3, 7), (3, 1, 3), (1, 2, 7)]
NAMES = [f'ax{ax}_p{pl}_c{c}' for (ax, pl, c) in SPEC]


def load_pool(L, include_existing=False):
    files = sorted(glob.glob(f'su3_k2pil_L{L}_N*_s*.npz'))
    if include_existing:
        f0 = f'rerun_su3_L{L}.npz'
        import os
        if os.path.exists(f0) and f0 not in files:
            files = [f0] + files
    if not files:
        return None, []
    d = {k: [] for k in ('a', 'H', 'G', 'P1', 'P2')}
    for fn in files:
        z = np.load(fn)
        for k in d:
            d[k].append(z[k])
    out = {k: np.concatenate(v, axis=0) for k, v in d.items()}
    return out, files


def mode_metrics(acc, j, rng, bw=0.30):
    a = acc['a'][:, j]
    kkde, skde = RC.kappa_kde(a, rng=rng)
    kap, seh, lam4, inv = RC.kappa_hist(a, rng=rng)
    k = RC.kernel_cond(a, {'H': acc['H'][:, j], 'G': acc['G'][:, j],
                           'P1': acc['P1'][:, j], 'P2': acc['P2'][:, j]}, bw=bw)
    EH, VG = k['E_H'], k['V_own']['G']
    EP2, VP1, cov = k['E_P2'], k['V_own']['P1'], k['C_G_P1']
    br = EH - VG
    W_geo = EP2 - VP1            # proven >= 0 (Phi convex => Prekopa)
    X = -2.0 * cov              # cross term (not sign-definite)
    W = W_geo + X               # total Jacobian contribution W_full
    # lag-1 autocorrelation of a (configs separated by `gap` sweeps)
    aa = a - a.mean()
    ac1 = float(np.sum(aa[1:] * aa[:-1]) / np.sum(aa * aa)) if a.size > 2 else np.nan
    return dict(gpp=kkde, se=skde, kh=kap, inv=inv, EH=EH, VG=VG, EP2=EP2,
                VP1=VP1, W_geo=W_geo, X=X, br=br, W=W, tot=br + W, ac1=ac1)


def report(L, include_existing):
    acc, files = load_pool(L, include_existing)
    lab = 'fresh+existing' if include_existing else 'fresh-only'
    if acc is None:
        print(f"SU(3) L={L} [{lab}]: no files yet")
        return None
    rng = np.random.default_rng(7000)
    N = acc['a'].shape[0]
    print(f"\n=== SU(3) L={L}  [{lab}]  N={N}  k=2pi/L={2*np.pi/L:.4f}  files={files}")
    print(f"  {'mode':14s} {'Gpp(KDE)':>9s} {'+-SE':>7s} {'kap_hist':>8s} {'1/Var':>6s} "
          f"{'E[H]':>6s} {'Var':>6s} {'E[Phi2]':>7s} {'Wgeo':>6s} {'X':>6s} {'W':>7s} {'br+W':>7s} {'ac1':>5s}")
    rows = []
    for j, nm in enumerate(NAMES):
        m = mode_metrics(acc, j, rng)
        rows.append(m)
        print(f"  {nm:14s} {m['gpp']:+9.4f} {m['se']:7.4f} {m['kh']:+8.4f} {m['inv']:6.3f} "
              f"{m['EH']:6.3f} {m['VG']:6.3f} {m['EP2']:7.3f} {m['W_geo']:+6.3f} "
              f"{m['X']:+6.3f} {m['W']:+7.3f} {m['tot']:+7.3f} {m['ac1']:5.2f}")
    gpp = np.array([m['gpp'] for m in rows]); se = np.array([m['se'] for m in rows])
    mean = gpp.mean(); se_tot = np.sqrt((se ** 2).mean() + gpp.var(ddof=1) / len(gpp))
    br = np.mean([m['br'] for m in rows]); W = np.mean([m['W'] for m in rows])
    Wgeo = np.mean([m['W_geo'] for m in rows]); X = np.mean([m['X'] for m in rows])
    tot = np.mean([m['tot'] for m in rows]); inv = np.mean([m['inv'] for m in rows])
    EP2 = np.mean([m['EP2'] for m in rows]); EH = np.mean([m['EH'] for m in rows])
    VG = np.mean([m['VG'] for m in rows]); VP1 = np.mean([m['VP1'] for m in rows])
    npos = int(np.sum(gpp > 0))
    print(f"  POOLED  Gpp={mean:+.3f}+-{se_tot:.3f}  ({npos}/6 modes>0; range "
          f"{gpp.min():+.3f}..{gpp.max():+.3f})  1/Var={inv:.3f}")
    print(f"          DECOMP E[H]={EH:.3f}  Var(G)={VG:.3f}  bracket={br:+.3f}  "
          f"E[Phi2]={EP2:.3f}  Var(Phi')={VP1:.3f}  Wgeo={Wgeo:+.3f}  X={X:+.3f}  "
          f"W={W:+.3f}  br+W={tot:+.3f}  (closure gpp-(br+W)={mean-tot:+.3f})")
    return dict(L=L, lab=lab, N=N, gpp=mean, se=se_tot, inv=inv, EH=EH, VG=VG,
                EP2=EP2, br=br, W=W, tot=tot)


def bw_sensitivity(L, include_existing):
    acc, _ = load_pool(L, include_existing)
    if acc is None:
        return
    lab = 'fresh+existing' if include_existing else 'fresh-only'
    print(f"\n  bandwidth sensitivity  SU(3) L={L} [{lab}]:")
    for bw in (0.20, 0.30, 0.40):
        rng = np.random.default_rng(7000)
        br = []; W = []; tot = []
        for j in range(acc['a'].shape[1]):
            m = mode_metrics(acc, j, rng, bw=bw)
            br.append(m['br']); W.append(m['W']); tot.append(m['tot'])
        print(f"    bw={bw:.2f}: bracket={np.mean(br):+.3f}  W={np.mean(W):+.3f}  "
              f"br+W={np.mean(tot):+.3f}")


if __name__ == '__main__':
    import os
    os.chdir('D:/Dev/math/research/unconditional/w7')
    out = []
    for L in (4, 6):
        r = report(L, include_existing=False)
        if r:
            out.append(r)
        r2 = report(L, include_existing=True)
        if r2:
            out.append(r2)
        bw_sensitivity(L, include_existing=False)
    print("\n=== SUMMARY (fresh+existing) ===")
    for i, L in enumerate((4, 6)):
        cand = [r for r in out if r['L'] == L and r['lab'] == 'fresh+existing']
        if cand:
            r = cand[0]
            print(f"  L={L}: N={r['N']}  Gpp={r['gpp']:+.3f}+-{r['se']:.3f}  "
                  f"[E[H]={r['EH']:.2f} Var={r['VG']:.2f} br={r['br']:+.3f} "
                  f"E[Phi2]={r['EP2']:.2f} W={r['W']:+.3f} br+W={r['tot']:+.3f}]")
