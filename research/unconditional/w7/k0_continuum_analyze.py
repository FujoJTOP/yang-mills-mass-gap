"""
K0-CONTINUUM analysis: G''(0) at k=2pi/L on ONE volume (L=8), with a
DELETE-BLOCK JACKKNIFE (robust to autocorrelation) for every quantity.

Exact decomposition (flat chart, conditioned at a=0):
    G''(0) = bracket + W
    bracket = E[H|0] - Var(G|0)              (action / dynamical part)
    W_geo   = E[Phi''|0] - Var(Phi'|0)       (Haar/Jacobian; PROVEN >=0 pointwise)
    X       = -2 Cov(G, Phi'|0)              (cross term)
    W       = W_geo + X
Gold-standard Gpp = curvature of -log P(a) (KDE), independent of the identity.

Output: per-mode table, per-k pooled table, overall pooled; jackknife SEs.
"""
import sys
import os
import json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rerun_core as RC                      # noqa: E402
from kzero_analyze import _kde_point         # noqa: E402


def stats(a, h, g, p1, p2, bw=0.30):
    """One sample's statistics (no bootstrap)."""
    k = RC.kernel_cond(a, {'H': h, 'G': g, 'P1': p1, 'P2': p2}, bw=bw)
    EH = k['E_H']; VG = k['V_own']['G']
    EP2 = k['E_P2']; VP1 = k['V_own']['P1']; cov = k['C_G_P1']
    return dict(Gpp=_kde_point(a), bracket=EH - VG,
                Wgeo=EP2 - VP1, X=-2.0 * cov,
                W=(EP2 - VP1) - 2.0 * cov,
                tot=(EH - VG) + (EP2 - VP1) - 2.0 * cov,
                EH=EH, VG=VG, EP2=EP2, VP1=VP1, cov=cov,
                inv=1.0 / np.var(a))


def delete_block_jack(a, h, g, p1, p2, nblk=20, bw=0.30):
    """Return (full_stat, se, mean_leave) for every key via delete-block
    jackknife over nblk contiguous blocks."""
    N = a.shape[0]
    nblk = int(min(nblk, max(4, N // 8)))
    edges = np.linspace(0, N, nblk + 1).astype(int)
    full = stats(a, h, g, p1, p2, bw)
    keys = list(full.keys())
    leave = {k: [] for k in keys}
    for b in range(nblk):
        keep = np.ones(N, bool)
        keep[edges[b]:edges[b + 1]] = False
        s = stats(a[keep], h[keep], g[keep], p1[keep], p2[keep], bw)
        for k in keys:
            leave[k].append(s[k])
    se = {}; ml = {}
    for k in keys:
        arr = np.array(leave[k])
        ml[k] = float(arr.mean())
        se[k] = float(np.sqrt((nblk - 1) / nblk * np.sum((arr - arr.mean()) ** 2)))
    return full, se, ml, nblk


def pooled_jack(A, H, G, P1, P2, nblk=20, bw=0.30):
    """Pooled over modes (equal weight) with delete-block jackknife."""
    N, Nm = A.shape
    nblk = int(min(nblk, max(4, N // 8)))
    edges = np.linspace(0, N, nblk + 1).astype(int)

    def pool(sl):
        rows = [stats(A[sl, j], H[sl, j], G[sl, j], P1[sl, j], P2[sl, j], bw)
                for j in range(Nm)]
        return {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}

    full = pool(slice(None))
    keys = list(full.keys())
    leave = {k: [] for k in keys}
    for b in range(nblk):
        keep = np.ones(N, bool)
        keep[edges[b]:edges[b + 1]] = False
        s = pool(keep)
        for k in keys:
            leave[k].append(s[k])
    se = {}; ml = {}
    for k in keys:
        arr = np.array(leave[k])
        ml[k] = float(arr.mean())
        se[k] = float(np.sqrt((nblk - 1) / nblk * np.sum((arr - arr.mean()) ** 2)))
    return full, se, ml, nblk


def fmt(name, full, se):
    return (f"  {name:20s} Gpp={full['Gpp']:+.4f}+-{se['Gpp']:.4f}  "
            f"br={full['bracket']:+.3f}+-{se['bracket']:.3f}  "
            f"Wgeo={full['Wgeo']:+.3f}+-{se['Wgeo']:.3f}  "
            f"X={full['X']:+.3f}+-{se['X']:.3f}  "
            f"W={full['W']:+.3f}+-{se['W']:.3f}  "
            f"br+W={full['tot']:+.3f}+-{se['tot']:.3f}  "
            f"1/Var={full['inv']:.3f}")


def main(path, nblk=20, bw=0.30, outdir=None):
    d = np.load(path)
    A, H, G, P1, P2 = d['a'], d['H'], d['G'], d['P1'], d['P2']
    N, Nm = A.shape
    names = [str(s) for s in d['names']] if 'names' in d else [f'm{j}' for j in range(Nm)]
    L = int(d['L']); beta = float(d['beta']); group = str(d['group'])
    kvals = np.asarray(d['kvals'], float) if 'kvals' in d else np.full(Nm, 2*np.pi/L)
    print(f"# {path}\n# {group} L={L} beta={beta} N={N} modes={Nm} "
          f"nblk={nblk} bw={bw}", flush=True)

    results = {'meta': dict(group=group, L=L, beta=beta, N=N, Nm=Nm,
                            nblk=nblk, bw=bw), 'per_mode': [], 'per_k': [], 'pooled': None}

    print("\n== per-mode (delete-block jackknife SE) ==", flush=True)
    for j in range(Nm):
        full, se, ml, nb = delete_block_jack(A[:, j], H[:, j], G[:, j], P1[:, j],
                                             P2[:, j], nblk=nblk, bw=bw)
        print(fmt(f"{names[j]} k={kvals[j]:.3f}", full, se), flush=True)
        results['per_mode'].append(dict(name=names[j], k=float(kvals[j]),
                                        full=full, se=se))

    print("\n== per-k pooled (equal weight over 6 modes at each n) ==", flush=True)
    uk = np.unique(np.round(kvals, 6))
    for kv in uk:
        js = [j for j in range(Nm) if abs(kvals[j] - kv) < 1e-6]
        full, se, ml, nb = pooled_jack(A[:, js], H[:, js], G[:, js],
                                       P1[:, js], P2[:, js], nblk=nblk, bw=bw)
        print(fmt(f"k={kv:.4f} ({len(js)} modes)", full, se), flush=True)
        results['per_k'].append(dict(k=float(kv), nmodes=len(js), full=full, se=se))

    print("\n== overall pooled (all modes) ==", flush=True)
    full, se, ml, nb = pooled_jack(A, H, G, P1, P2, nblk=nblk, bw=bw)
    print(fmt("ALL", full, se), flush=True)
    # sign / dominance summary
    frac = full['W'] / full['Gpp'] if full['Gpp'] != 0 else float('nan')
    print(f"\n  SIGN: Gpp={full['Gpp']:+.4f}+-{se['Gpp']:.4f}  "
          f"({full['Gpp']/se['Gpp']:+.1f} sigma) ; "
          f"W carries {100*frac:.0f}% of Gpp ; "
          f"Wgeo={full['Wgeo']:+.3f}+-{se['Wgeo']:.3f} ; "
          f"closure br+W-Gpp={full['tot']-full['Gpp']:+.3f}", flush=True)
    results['pooled'] = dict(full=full, se=se, W_frac_of_Gpp=float(frac))

    if outdir is None:
        outdir = HERE
    stem = os.path.splitext(os.path.basename(path))[0]
    jout = os.path.join(outdir, f'{stem}_results.json')
    with open(jout, 'w') as f:
        json.dump(results, f, indent=1)
    print(f"\n  wrote {jout}", flush=True)
    return results


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('file')
    ap.add_argument('--nblk', type=int, default=20)
    ap.add_argument('--bw', type=float, default=0.30)
    ap.add_argument('--outdir', default=None)
    a = ap.parse_args()
    main(a.file, nblk=a.nblk, bw=a.bw, outdir=a.outdir)
