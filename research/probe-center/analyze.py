"""Analyse probe_*.npz: vortex structure, signed/absolute ratios, sector tests.

Prints tables and writes plots to plots/.
"""
import glob
import os
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = 'D:/Dev/math/research/probe-center'
os.makedirs(f'{HERE}/plots', exist_ok=True)


def load_all():
    out = {}
    for fn in sorted(glob.glob(f'{HERE}/probe_b*.npz')):
        lbl = re.sub(r'probe_|\.npz', '', os.path.basename(fn))
        d = np.load(fn, allow_pickle=True)
        if d['nsnap'] == 0:
            continue
        out[lbl] = d
    return out


def jack_err(x, axis=0):
    x = np.asarray(x, float)
    n = x.shape[axis]
    if n < 2:
        return np.zeros_like(np.nanmean(x, axis=axis))
    m = np.nanmean(x, axis=axis)
    s = np.nanstd(x, axis=axis, ddof=1) / np.sqrt(n)
    return s


def creutz(W):
    """Creutz ratio chi(R,T) = -ln[W(R,T)W(R-1,T-1)/(W(R,T-1)W(R-1,T))]."""
    n = W.shape[0]
    C = np.full((n, n), np.nan)
    for R in range(2, n):
        for T in range(2, n):
            num = W[R, T] * W[R - 1, T - 1]
            den = W[R, T - 1] * W[R - 1, T]
            if num > 0 and den > 0:
                C[R, T] = -np.log(num / den)
    return C


def main():
    data = load_all()
    print("=" * 78)
    print("VORTEX STRUCTURE  (mean over snapshots +- s.e.m.)")
    print("=" * 78)
    hdr = (f"{'label':10s} {'L':>2s} {'beta':>5s} {'<R>':>7s} {'dens':>7s} "
           f"{'nvort':>7s} {'ncl':>5s} {'k1':>6s} {'k2':>6s} {'maxsz':>7s}")
    print(hdr)
    rows = {}
    for lbl, d in data.items():
        L = int(d['L']); beta = float(d['beta'])
        dens = d['density']; nv = d['nvortex']
        k1 = d['k1']; k2 = d['k2']; ncl = d['nclust']
        allsz = d['allsizes']
        maxsz = allsz.max() if allsz.size else 0
        rows[lbl] = dict(L=L, beta=beta,
                         R=(d['R'].mean(), jack_err(d['R'])),
                         dens=(dens.mean(), jack_err(dens)),
                         nv=(nv.mean(), jack_err(nv)),
                         ncl=(ncl.mean(), jack_err(ncl)),
                         k1=(k1.mean(), jack_err(k1)),
                         k2=(k2.mean(), jack_err(k2)),
                         maxsz=maxsz, allsz=allsz,
                         d=d)
        print(f"{lbl:10s} {L:2d} {beta:5.3f} {d['R'].mean():7.3f} "
              f"{dens.mean():7.3f} {nv.mean():7.1f} {ncl.mean():5.1f} "
              f"{k1.mean():6.1f} {k2.mean():6.1f} {maxsz:7d}")

    print()
    print("=" * 78)
    print("SECTOR TEST 1: fundamental Wilson loop, full vs center-projected")
    print("  W_F = Re<Tr/3>;  adjoint W_A = (9<|Tr/3|^2>-1)/8")
    print("=" * 78)
    for lbl, r in rows.items():
        d = r['d']; Tm = d['Tf'].mean(0); T2 = d['Tf2'].mean(0)
        Tp = d['Tp'].mean(0); Tp2 = d['Tp2'].mean(0)
        WF = Tm.real
        WA = (9 * T2 - 1) / 8
        WFp = Tp.real
        WAp = (9 * Tp2 - 1) / 8
        print(f"\n--- {lbl}  L={r['L']} beta={r['beta']} ---")
        print("  R,T |   W_F full |  W_F proj |  W_A full |  W_A proj")
        n = WF.shape[0]
        for R in range(1, n):
            for T in range(1, n):
                print(f"  {R},{T} | {WF[R,T]:10.4f} | {WFp[R,T]:9.4f} | "
                      f"{WA[R,T]:9.4f} | {WAp[R,T]:9.5f}")
        # Creutz ratios on the fundamental
        cF = creutz(WF); cFp = creutz(WFp)
        print("  Creutz fundamental full:",
              np.array2string(cF, precision=3, suppress_small=True))
        print("  Creutz fundamental proj:",
              np.array2string(cFp, precision=3, suppress_small=True))
        print("  adjoint W_A(full) 1,1 / 2,2:", f"{WA[1,1]:.4f} / {WA[2,2]:.4f}")
        print("  adjoint W_A(proj) max deviation from 1:",
              f"{np.abs(WAp[1:,1:]-1).max():.2e}")

    print()
    print("=" * 78)
    print("SECTOR TEST 2: Polyakov correlators (charged = fundamental, neutral = adjoint)")
    print("=" * 78)
    for lbl, r in rows.items():
        d = r['d']
        Cf = d['Cf'].mean(0) / 9.0      # fundamental, normalised
        CA = d['CA'].mean(0)            # adjoint connected (real)
        Cfp = d['Cfp'].mean(0) / 9.0
        print(f"\n--- {lbl} ---")
        print("  r | C_F full(Re) | C_F proj(Re) | C_A full(conn) | <|P|>full  <|P|>proj")
        nmax = len(Cf)
        for rr in range(nmax):
            print(f"  {rr} | {Cf[rr].real:12.4e} | {Cfp[rr].real:12.4e} | "
                  f"{CA[rr]:13.4e} | {d['poly_mean_full'].mean():.4f}   "
                  f"{d['poly_mean_proj'].mean():.4f}")

    print()
    print("=" * 78)
    print("SIGNED / ABSOLUTE ACTIVITY RATIO  z_signed/z_abs")
    print("  D1: sign = Re(omega^k)/|omega^k|  (physical real phase; k in {1,2})")
    print("  D2: sign = orientation  (k=1:+1, k=2:-1)")
    print("=" * 78)
    for lbl, r in rows.items():
        d = r['d']
        k1 = r['k1'][0]; k2 = r['k2'][0]
        # D1: every nontrivial plaquette has Re(omega^k) = -1/2 -> ratio = -0.5
        D1 = -0.5
        # D2 over the whole ensemble
        D2 = (k1 - k2) / max(k1 + k2, 1)
        print(f"{lbl:10s}: D1 = {D1:+.4f} (size-independent, group theory)   "
              f"D2(ensemble) = {D2:+.4f}   nplaquettes={int(k1+k2)}")
    # size-binned D2 from matched per-cluster (size, N1-N2) arrays
    print("\n  Size-binned z_signed/z_abs (orientation convention, per cluster):")
    print("   bin(s)        | N_clusters | <(N1-N2)/(N1+N2)> | asym. |N1-N2|/(N1+N2)")
    bins = [(1, 1), (2, 3), (4, 7), (8, 15), (16, 63), (64, 255), (256, 10**9)]
    for lbl, r in rows.items():
        d = r['d']
        if 'cl_size' not in d.files:
            continue
        sz = d['cl_size'].astype(float); cq = d['cq_all'].astype(float)
        out = []
        for lo, hi in bins:
            m = (sz >= lo) & (sz <= hi)
            if m.sum() == 0:
                out.append(f"{lo}-{hi}:--")
                continue
            ratio = cq[m] / sz[m]
            out.append(f"{lo}-{hi}: n={int(m.sum())} <r>={ratio.mean():+.3f} "
                       f"<|r|>={np.abs(ratio).mean():.3f}")
        print(f"  {lbl}:")
        for o in out:
            print("     ", o)

    # ---------- plots ----------
    # 1) density vs beta
    fig, ax = plt.subplots(figsize=(6, 4))
    for L in sorted({r['L'] for r in rows.values()}):
        xs = [r['beta'] for r in rows.values() if r['L'] == L]
        ys = [r['dens'][0] for r in rows.values() if r['L'] == L]
        ye = [r['dens'][1] for r in rows.values() if r['L'] == L]
        o = np.argsort(xs)
        ax.errorbar(np.array(xs)[o], np.array(ys)[o], yerr=np.array(ye)[o],
                    marker='o', label=f'L={L}')
    ax.set_xlabel('beta'); ax.set_ylabel('vortex density / plaquette')
    ax.set_title('Z(3)-projected vortex density (MCG)')
    ax.legend(); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(f'{HERE}/plots/vortex_density.png', dpi=130)
    plt.close(fig)

    # 2) W_F full vs proj and W_A for a representative config
    for lbl in ['b2.6_L6', 'b3.0_L6', 'b6.0_L6']:
        if lbl not in rows:
            continue
        d = rows[lbl]['d']
        Tm = d['Tf'].mean(0); Tp = d['Tp'].mean(0); T2 = d['Tf2'].mean(0)
        WF = Tm.real; WFp = Tp.real; WA = (9*T2-1)/8
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        ax[0].plot(range(1, WF.shape[1]), WF[1, 1:], 'o-', label='full fund.')
        ax[0].plot(range(1, WFp.shape[1]), WFp[1, 1:], 's--', label='proj fund.')
        ax[0].plot(range(1, WA.shape[1]), WA[1, 1:], '^-', label='full adjoint')
        ax[0].set_xlabel('T (R=1)'); ax[0].set_ylabel('W')
        ax[0].set_yscale('symlog', linthresh=1e-3); ax[0].legend(); ax[0].grid(alpha=.3)
        ax[0].set_title(f'{lbl}: W(R=1,T)')
        ax[1].plot(range(1, WF.shape[0]), WF[1:, 1], 'o-', label='full fund.')
        ax[1].plot(range(1, WFp.shape[0]), WFp[1:, 1], 's--', label='proj fund.')
        ax[1].plot(range(1, WA.shape[0]), WA[1:, 1], '^-', label='full adjoint')
        ax[1].set_xlabel('R (T=1)'); ax[1].set_ylabel('W')
        ax[1].set_yscale('symlog', linthresh=1e-3); ax[1].legend(); ax[1].grid(alpha=.3)
        ax[1].set_title(f'{lbl}: W(R,T=1)')
        fig.tight_layout(); fig.savefig(f'{HERE}/plots/wilson_{lbl}.png', dpi=130)
        plt.close(fig)

    # 3) Polyakov correlators
    for lbl in ['b2.6_L6', 'b3.0_L6', 'b6.0_L6']:
        if lbl not in rows:
            continue
        d = rows[lbl]['d']
        Cf = d['Cf'].mean(0).real / 9.0
        Cfp = d['Cfp'].mean(0).real / 9.0
        CA = d['CA'].mean(0)
        fig, ax = plt.subplots(figsize=(6, 4))
        r = np.arange(len(Cf))
        ax.plot(r, Cf, 'o-', label='fund. full')
        ax.plot(r, Cfp, 's--', label='fund. proj')
        ax.plot(r, CA / np.abs(CA[0]) if CA[0] != 0 else CA, '^-', label='adjoint conn. (norm.)')
        ax.set_xlabel('r'); ax.set_ylabel('C(r)')
        ax.set_yscale('symlog', linthresh=1e-5); ax.legend(); ax.grid(alpha=.3)
        ax.set_title(f'{lbl}: Polyakov correlators')
        fig.tight_layout(); fig.savefig(f'{HERE}/plots/polyakov_{lbl}.png', dpi=130)
        plt.close(fig)

    # 4) cluster size histogram (finite clusters) across configs
    fig, ax = plt.subplots(figsize=(6, 4))
    for lbl, r in rows.items():
        s = r['allsz']
        if s.size == 0:
            continue
        hist, edges = np.histogram(s, bins=np.logspace(0, np.log10(max(2, s.max()) + 1), 15))
        centers = 0.5 * (edges[1:] + edges[:-1])
        m = hist > 0
        ax.loglog(centers[m], hist[m], 'o-', label=lbl)
    ax.set_xlabel('cluster size (plaquettes)'); ax.set_ylabel('count')
    ax.set_title('vortex cluster size distribution')
    ax.legend(fontsize=7); ax.grid(alpha=.3, which='both')
    fig.tight_layout(); fig.savefig(f'{HERE}/plots/cluster_sizes.png', dpi=130)
    plt.close(fig)

    print("\nplots written to", f'{HERE}/plots')


if __name__ == '__main__':
    main()
