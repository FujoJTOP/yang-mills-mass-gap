"""Summarise kappa_results_*.npz into a CSV + plot."""
import glob
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    recs = []
    for f in sorted(glob.glob(os.path.join(HERE, 'kappa_results_*.npz'))):
        d = dict(np.load(f, allow_pickle=True))
        L = int(d['L'][0]) if d['L'].size else None
        betas = np.unique(d['beta'])
        for beta in betas:
            k = d['kappa'][d['beta'] == beta]
            r = d['lanczos_res'][d['beta'] == beta]
            hmin = d['hl_min'][d['beta'] == beta]
            hneg = d['hl_neg_frac'][d['beta'] == beta]
            hmean = d['hl_mean'][d['beta'] == beta]
            recs.append(dict(L=L, beta=beta, kappa_mean=float(k.mean()),
                             kappa_std=float(k.std(ddof=1)) if k.size > 1 else 0.0,
                             kappa_min=float(k.min()), kappa_max=float(k.max()),
                             frac_negative=float(np.mean(k < 0)),
                             lanczos_res_max=float(r.max()),
                             hl_min=float(hmin.min()),
                             hl_neg_frac=float(hneg.mean()),
                             hl_mean=float(hmean.mean()), n=int(k.size)))

    with open(os.path.join(HERE, 'kappa_summary.csv'), 'w') as f:
        f.write('L,beta,kappa_mean,kappa_std,kappa_min,kappa_max,frac_negative,'
                'lanczos_res_max,hl_min,hl_neg_frac,hl_mean,n_snapshots\n')
        for r in recs:
            f.write(f"{r['L']},{r['beta']:.4f},{r['kappa_mean']:.6f},"
                    f"{r['kappa_std']:.6f},{r['kappa_min']:.6f},{r['kappa_max']:.6f},"
                    f"{r['frac_negative']:.3f},{r['lanczos_res_max']:.3g},"
                    f"{r['hl_min']:.5f},{r['hl_neg_frac']:.4f},{r['hl_mean']:.5f},"
                    f"{r['n']}\n")
    print(open(os.path.join(HERE, 'kappa_summary.csv')).read())

    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for L in sorted({r['L'] for r in recs}):
        rr = sorted([r for r in recs if r['L'] == L], key=lambda x: x['beta'])
        xs = [r['beta'] for r in rr]
        axs[0].errorbar(xs, [r['kappa_mean'] for r in rr],
                        yerr=[r['kappa_std'] for r in rr], marker='o', capsize=3,
                        label=f'L={L}')
        axs[1].plot(xs, [r['hl_neg_frac'] for r in rr], 'o-', label=f'L={L}')
    axs[0].axhline(0, color='k', lw=.8)
    axs[0].set_xlabel(r'$\beta$'); axs[0].set_ylabel(r'$\kappa_{\rm eff}$')
    axs[0].set_title('lowest quotient-Hessian eigenvalue'); axs[0].legend()
    axs[0].grid(alpha=.3)
    axs[1].set_xlabel(r'$\beta$'); axs[1].set_ylabel('fraction of links with $\\lambda_{\\min}<0$')
    axs[1].set_title('local 8x8 link curvature'); axs[1].legend(); axs[1].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, 'plots', 'kappa_eff.png'), dpi=140)


if __name__ == '__main__':
    main()
