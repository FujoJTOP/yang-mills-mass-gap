"""Distribution figure: true blocked ensemble vs fitted model vs Wilson-only."""
import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import imc
import blk


def hist_of(data, fit, sweeps=1500):
    d = dict(np.load(data))
    f = json.load(open(fit))
    L = int(d['L']) // 2
    h_true = np.sum(np.array(d['hist']), axis=0)
    recs, a, _ = imc.model_run(np.array(f['lam']), f['ops'], L, sweeps, 400, 2,
                               11, verbose=False, full=True)
    h_mod = np.sum(np.array([r['hist'] for r in recs]), axis=0)
    return d, f, h_true, h_mod


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--fit', required=True)
    ap.add_argument('--ponly', type=float, default=0.0,
                    help='also show a Wilson-only model at this lambda_P')
    ap.add_argument('--out', default='plots/blocked_dist.pdf')
    args = ap.parse_args()
    os.makedirs('plots', exist_ok=True)
    d, f, h_true, h_mod = hist_of(args.data, args.fit)
    edges = blk.HIST_EDGES
    ctr = 0.5 * (edges[:-1] + edges[1:])
    w = edges[1] - edges[0]
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.4))
    for a in ax:
        a.step(ctr, h_true / h_true.sum() / w, where='mid', lw=1.6,
               label='blocked ensemble')
        a.step(ctr, h_mod / h_mod.sum() / w, where='mid', lw=1.2,
               label=f"fit {f['ops']} $\\lambda$={np.round(f['lam'],3)}")
    if args.ponly > 0:
        L = int(d['L']) // 2
        recs, _, _ = imc.model_run([args.ponly], ['P'], L, 1500, 400, 2, 13,
                                   verbose=False, full=True)
        h_p = np.sum(np.array([r['hist'] for r in recs]), axis=0)
        for a in ax:
            a.step(ctr, h_p / h_p.sum() / w, where='mid', lw=1.0, ls='--',
                   label=f'Wilson only $\\lambda_P$={args.ponly}')
    ax[0].set_xlim(-0.35, 0.55)
    ax[0].set_xlabel(r'blocked plaquette  $\frac13$Re tr $V_p$')
    ax[0].set_ylabel('density')
    ax[0].legend(fontsize=6)
    ax[1].semilogy(ctr, h_true / h_true.sum() / w, '.', ms=2)
    ax[1].semilogy(ctr, h_mod / h_mod.sum() / w, '.', ms=2)
    if args.ponly > 0:
        ax[1].semilogy(ctr, h_p / h_p.sum() / w, '.', ms=2)
    ax[1].set_xlim(-0.4, 0.4)
    ax[1].set_xlabel(r'$\frac13$Re tr $V_p$')
    ax[1].set_ylabel('density (log)')
    fig.tight_layout()
    fig.savefig(args.out)
    print('wrote', args.out)


if __name__ == '__main__':
    main()
