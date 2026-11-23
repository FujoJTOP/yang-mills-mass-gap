"""Final SU(3) beta=5.7 k=2pi/L analysis with autocorrelation-robust block bootstrap.

Pools (a,H,G,P1,P2) per L from:
  - fresh  su3_k2pil_L{L}_N*_s*.npz
  - optionally the pre-existing rerun_su3_L{L}.npz
injects L/beta/kphys/names metadata, and runs kzero_analyze.analyse (block
bootstrap, block=20, 400 resamples) to get G''(0)=curvature of -log P(a) and the
exact-identity decomposition  G''=(E[H]-Var(G))+W  at a=0 (bw=0.30).
"""
import sys
import glob
import os
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import kzero_analyze as KZ     # noqa: E402

SPEC = [(1, 2, 3), (1, 3, 3), (0, 2, 3), (2, 3, 7), (3, 1, 3), (1, 2, 7)]
NAMES = [f'ax{ax}_p{pl}_c{c}' for (ax, pl, c) in SPEC]
KEYS = ('a', 'H', 'G', 'P1', 'P2')


def build_pool(L, beta=5.7, include_existing=False, out=None):
    files = sorted(glob.glob(f'su3_k2pil_L{L}_N*_s*.npz'))
    if include_existing:
        f0 = f'rerun_su3_L{L}.npz'
        if os.path.exists(f0):
            files = [f0] + files
    if not files:
        return None, []
    acc = {k: [] for k in KEYS}
    for fn in files:
        z = np.load(fn)
        for k in KEYS:
            acc[k].append(z[k])
    acc = {k: np.concatenate(v, axis=0) for k, v in acc.items()}
    meta = dict(L=L, beta=beta, kphys=2 * np.pi / L,
                names=np.array(NAMES), files=np.array(files))
    if out:
        np.savez(out, **acc, **meta)
    return (acc, meta), files


if __name__ == '__main__':
    os.chdir('D:/Dev/math/research/unconditional/w7')
    mode = sys.argv[1] if len(sys.argv) > 1 else 'fresh'
    inc = (mode == 'pool')
    print("#" * 90)
    print(f"# SU(3) beta=5.7, softest nonzero mode k=2pi/L, {'fresh+existing' if inc else 'fresh-only'}")
    print("#" * 90)
    summary = []
    for L in (4, 6):
        built, files = build_pool(L, include_existing=inc, out=f'_pool_su3_L{L}.npz')
        if built is None:
            print(f"\nSU(3) L={L}: NO FILES (fresh run not finished)\n")
            continue
        acc, meta = built
        N = acc['a'].shape[0]
        print(f"\n>>> SU(3) L={L}  N={N}  k=2pi/L={2*np.pi/L:.4f}  files={[os.path.basename(f) for f in files]}")
        rows, pooled = KZ.analyse(f'_pool_su3_L{L}.npz', bw=0.30, nboot=400, block=20,
                                  rngseed=7000, verbose=True)
        tau = np.array([r['tau'] for r in rows])
        npos = int(sum(1 for r in rows if r['kap'] > 0))
        # CONSERVATIVE pooled SE: the 6 modes share configs (correlated), so do
        # NOT divide the bootstrap term by sqrt(Nm); keep one-mode SE scale.
        se_single = float(np.sqrt(np.mean([r['se'] ** 2 for r in rows])))
        se_mode = float(np.std([r['kap'] for r in rows], ddof=1) / np.sqrt(len(rows)))
        se_cons = float(np.sqrt(se_single ** 2 + se_mode ** 2))
        print(f"  >>> SIGN: {npos}/{len(rows)} modes have G''>0 ; pooled G''={pooled['kap']:+.4f}"
              f" +- {se_cons:.4f} (conservative; one-mode block-boot scale) "
              f"[kzero √Nm SE={pooled['se']:.4f}]  (mean tau={tau.mean():.2f} stored-config lags)")
        summary.append(dict(L=L, N=N, gpp=pooled['kap'], se=se_cons, se_kz=pooled['se'],
                            br=pooled['bracket'], W=pooled['W'], tot=pooled['tot'],
                            EH=pooled['EH'], VG=pooled['VG'], EP2=pooled['EP2'],
                            VP1=pooled['VP1'], cov=pooled['cov'], inv=pooled['inv'],
                            npos=npos, nm=len(rows), tau=tau.mean()))
    print("\n" + "=" * 90)
    print(f"FINAL SUMMARY [{mode}]")
    print(f"{'L':>2s} {'N':>5s} {'Gpp':>8s} {'+-SE':>7s} {'E[H]':>6s} {'Var':>6s} "
          f"{'bracket':>8s} {'E[Phi2]':>7s} {'W_geo':>6s} {'X':>6s} {'W':>7s} {'br+W':>7s} {'sign':>6s}")
    for s in summary:
        wgeo = s['EP2'] - s['VP1']
        X = -2.0 * s['cov']
        print(f"{s['L']:>2d} {s['N']:>5d} {s['gpp']:+8.3f} {s['se']:7.3f} {s['EH']:6.3f} "
              f"{s['VG']:6.3f} {s['br']:+8.3f} {s['EP2']:7.3f} {wgeo:+6.3f} {X:+6.3f} "
              f"{s['W']:+7.3f} {s['tot']:+7.3f} {s['npos']:>2d}/{s['nm']:<2d}")
    print("=" * 90)
