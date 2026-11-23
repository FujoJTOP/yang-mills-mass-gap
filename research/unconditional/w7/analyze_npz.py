"""Re-analyse saved per-mode data (partA_*_N*.npz): kappa_hist, bootstrap SE,
quartic, 1/Var, bare E[H]/Var(G), and conditional width-scan Ev(w) that should
approach kappa_hist as the window narrows (identity + binning check)."""
import sys
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import gpp_core as C   # noqa: E402
from w7_final import width_scan   # noqa: E402


def modes_of(d):
    names = []
    for k in d.files:
        if k.endswith('_a'):
            names.append(k[:-2])
    return names


def main(fn, safe_key=None, safe_thresh=None):
    d = np.load(fn)
    names = modes_of(d)
    ang = d['ang'] if 'ang' in d.files else None
    safe = None
    if ang is not None and safe_thresh is not None:
        safe = ang < safe_thresh
        print(f"# {fn}: {safe.sum()}/{ang.size} branch-safe (max|ang|<{safe_thresh:.3f})")
    else:
        print(f"# {fn}: {len(names)} modes, N={len(d[names[0]+'_a'])}")
    rng = np.random.default_rng(12345)
    hdr = (f"{'mode':22s} {'Gpp':>8s} {'+-SE':>7s} {'1/Var':>7s} {'E[H]':>7s} "
           f"{'VarG':>7s} {'lam4':>7s} {'sd':>6s}   width-scan Ev(w)=E[H|w]-Var(G|w)")
    print(hdr)
    for n in names:
        a = d[n + '_a']; H = d[n + '_H']; G = d[n + '_G']
        kap, se, lam, _, inv = C.veff_curvature(a, rng=rng, nboot=500)
        sc = width_scan(a, H, G)
        scs = " ".join(f"{w:.2f}:{eh-vg:+.3f}" if not np.isnan(eh) else f"{w}:nan"
                       for (w, eh, vg, cnt) in sc)
        extra = ""
        if safe is not None:
            ks, ses, *_ = C.veff_curvature(a[safe], rng=rng, nboot=400)
            extra = f"  |safe:{ks:+.3f}+-{ses:.3f}"
        print(f"{n:22s} {kap:+8.4f} {se:7.4f} {inv:+7.4f} {H.mean():7.3f} "
              f"{G.var():7.3f} {lam:+7.3f} {a.std():6.3f}   {scs}{extra}")


if __name__ == '__main__':
    fn = sys.argv[1]
    sk = 'ang' if (len(sys.argv) > 2) else None
    if sk:
        main(fn, safe_key=sk, safe_thresh=float(sys.argv[2]) if len(sys.argv) > 2 else 2.99)
    else:
        main(fn)
