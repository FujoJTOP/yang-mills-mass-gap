"""Fresh SU(3) beta=5.7 runs at the SOFTEST nonzero mode k=2pi/L (n=1).

Reuses rerun_run.measure / analyse_mode (flat exponential chart, branch-safe
traceless log, analytic Jacobian derivatives).  Saves to a NEW filename so the
existing rerun_su3_L*.npz are never overwritten.

For each L: 6 transverse plane-wave modes at n=1  =>  |k|=2pi/L.
Reports, per mode and pooled:
  Gpp(KDE)  = curvature of -log P(a)            (gold standard, exact identity RHS)
  bracket   = E[H]-Var(G)  at a=0               (action / dynamical piece)
  E[Phi'']  , W = E[Phi'']-Var(Phi')-2Cov(G,Phi')
  bracket+W = identity RHS, must close to Gpp(KDE)
  1/Var(a)
"""
import sys
import time
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import rerun_run as R          # noqa: E402
import gpp_core3 as G3         # noqa: E402
import su3 as S3               # noqa: E402

SPEC = [(1, 2, 3), (1, 3, 3), (0, 2, 3), (2, 3, 7), (3, 1, 3), (1, 2, 7)]
NAMES = [f'ax{ax}_p{pl}_c{c}' for (ax, pl, c) in SPEC]


def main(L, N, seed, therm, gap):
    beta = 5.7
    print(f"# FRESH SU(3) L={L} beta={beta} N={N} seed={seed} therm={therm} gap={gap} "
          f"k=2pi/L={2*np.pi/L:.4f}", flush=True)
    t0 = time.time()
    cfgs = R.sample_su3(L, beta, therm, N, gap, seed)
    tsamp = time.time() - t0
    plaq = np.mean([S3.mean_plaquette(U) for U in cfgs[:min(200, len(cfgs))]])
    print(f"  sampled {len(cfgs)} cfgs in {tsamp:.0f}s; mean plaquette(200)={plaq:.4f}",
          flush=True)
    modes = [G3.plane_mode(L, ax, pl, c, n=1, kind='cos') for (ax, pl, c) in SPEC]
    acc = {k: [] for k in ('a', 'H', 'G', 'P1', 'P2')}
    tw = []
    t0 = time.time()
    for i, U in enumerate(cfgs):
        r = R.measure(U, beta, modes, 'su3')
        for k in acc:
            acc[k].append(r[k])
        tw.append(r['ntwist'])
        if (i + 1) % 250 == 0:
            print(f"    measured {i+1}/{N} ({time.time()-t0:.0f}s)", flush=True)
    for k in acc:
        acc[k] = np.array(acc[k])
    tw = np.array(tw)
    nlink = 4 * L ** 4
    print(f"  center-twist: {tw.sum()}/{len(tw)*nlink} links "
          f"({100*tw.sum()/(len(tw)*nlink):.2f}%); corrected log exact for ALL", flush=True)
    out = f'su3_k2pil_L{L}_N{N}_s{seed}.npz'
    np.savez(out, **acc)
    print(f"  saved {out}", flush=True)
    return acc


def analyse(acc, L, label):
    rng = np.random.default_rng(7000)
    A = acc['a']
    Nm = A.shape[1]
    rows = []
    print(f"\n  {label}: N={A.shape[0]}  modes={Nm}  k=2pi/L={2*np.pi/L:.4f}")
    hdr = (f"  {'mode':14s} {'Gpp(KDE)':>9s} {'+-SE':>7s} {'kap_hist':>8s} {'1/Var':>6s} "
           f"{'E[H]':>6s} {'VarG':>6s} {'E[Phi2]':>7s} {'bracket':>8s} {'W':>7s} {'br+W':>7s}")
    print(hdr, flush=True)
    for j, nm in enumerate(NAMES):
        r = R.analyse_mode(A[:, j], acc['H'][:, j], acc['G'][:, j],
                           acc['P1'][:, j], acc['P2'][:, j], rng, bw=0.30)
        c = r['cond']
        rows.append(dict(name=nm, gpp=r['kkde'], se=r['skde'], kh=r['kap'],
                         inv=r['inv'], EH=c['EH'], VG=c['VG'], EP2=c['EP2'],
                         br=c['bracket'], W=c['W'], tot=c['tot']))
        print(f"  {nm:14s} {r['kkde']:+9.4f} {r['skde']:7.4f} {r['kap']:+8.4f} "
              f"{r['inv']:6.3f} {c['EH']:6.3f} {c['VG']:6.3f} {c['EP2']:7.3f} "
              f"{c['bracket']:+8.3f} {c['W']:+7.3f} {c['tot']:+7.3f}", flush=True)
    kk = np.array([x['gpp'] for x in rows]); se = np.array([x['se'] for x in rows])
    mean = kk.mean()
    se_tot = np.sqrt((se ** 2).mean() + kk.var(ddof=1) / len(kk))
    br = np.mean([x['br'] for x in rows]); W = np.mean([x['W'] for x in rows])
    tot = np.mean([x['tot'] for x in rows]); inv = np.mean([x['inv'] for x in rows])
    EP2 = np.mean([x['EP2'] for x in rows]); EH = np.mean([x['EH'] for x in rows])
    VG = np.mean([x['VG'] for x in rows])
    print(f"  POOLED {label:10s} Gpp={mean:+.3f}+-{se_tot:.3f} (range {kk.min():+.3f}..{kk.max():+.3f}) "
          f"1/Var={inv:.3f} | E[H]={EH:.3f} Var={VG:.3f} bracket={br:+.3f} "
          f"E[Phi2]={EP2:.3f} W={W:+.3f} br+W={tot:+.3f}", flush=True)
    return dict(L=L, N=A.shape[0], gpp=mean, se=se_tot, inv=inv, EH=EH, VG=VG,
                EP2=EP2, br=br, W=W, tot=tot)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, required=True)
    ap.add_argument('--N', type=int, required=True)
    ap.add_argument('--seed', type=int, default=303)
    ap.add_argument('--therm', type=int, default=600)
    ap.add_argument('--gap', type=int, default=0)
    a = ap.parse_args()
    gap = a.gap if a.gap > 0 else (4 if a.L <= 4 else 5)
    acc = main(a.L, a.N, a.seed, a.therm, gap)
    analyse(acc, a.L, f'SU(3) L={a.L}')
