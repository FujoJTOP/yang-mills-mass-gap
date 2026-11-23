"""NA-3-T analysis -> results_eta.csv + na3t_table.txt.

Every grid point is labelled `cert` or `est`:
  cert = a recorded witness pair (config U, link e, 18 neighbour replacements g)
         yields a lower bound  sum_e' TV_lb >= 1  on the Dobrushin row sum.
         Since  eta >= sum_e' I(e,e') >= sum_e' TV_lb,  this is a strict kill.
  est  = no witness cleared 1; the point estimate is only evidence, never a
         survival certificate.
TV_lb is the witness-function bound  (1/2)|E_p f - E_q f|  with
f(U)=Re tr(U dS)/||dS||_* in [-1,1], minus 4 sigma of Monte-Carlo error.
"""
import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TAU_L4 = {2.2: 0.52, 3.0: 0.50, 6.0: 1.61}


def main():
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, 'runs', 'eta_*.npz'))):
        d = dict(np.load(f, allow_pickle=True))
        L = int(d['L']); beta = float(d['beta'])
        er = np.asarray(d['eta_rows'], float)
        lb = np.asarray(d.get('lb_rows', np.zeros_like(er)), float)
        ub = np.asarray(d['ub_rows'], float)
        n = len(er)
        pairs = np.asarray(d['tv_pairs'], float)
        lb_max = float(lb.max()) if n else np.nan
        cert = bool(lb_max >= 1.0)
        wname = f'witness_L{L}_b{beta:.1f}.npz'
        rows.append(dict(L=L, beta=beta, op='link', eta=float(er.mean()),
                         eta_err=float(er.std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0,
                         eta_max=float(er.max()), ub=float(ub.mean()),
                         lb_max=lb_max, lb_mean=float(lb.mean()),
                         tv_mean=float(pairs[:, 1].mean()) if len(pairs) else np.nan,
                         ds_mean=float(pairs[:, 0].mean()) if len(pairs) else np.nan,
                         n_cfg=n, nsamp=int(d['nsamp']),
                         cert='cert' if cert else 'est',
                         wname=wname if cert else ''))

    with open(os.path.join(HERE, 'results_eta.csv'), 'w') as f:
        f.write('L,beta,op,eta,eta_err,lam,lam_err,tau_hb,m_eff,ratio,n_cfg,'
                'kernel,cert,notes\n')
        for r in rows:
            notes = (f'eta=typical-config row sum of BAR-TV; witness lb '
                     f'{r["lb_max"]:.2f} (mean {r["lb_mean"]:.2f}); '
                     f'rigorous log-density ub sum={r["ub"]:.2f}')
            if r['cert']:
                notes += f'; WITNESS={r["wname"]}'
            f.write(f"{r['L']},{r['beta']:.3f},{r['op']},{r['eta']:.4f},"
                    f"{r['eta_err']:.4f},nan,nan,"
                    f"{TAU_L4.get(round(r['beta'], 1), float('nan')):.2f},nan,nan,"
                    f"{r['n_cfg']},link_CM,{r['cert']},\"{notes}\"\n")
    print(open(os.path.join(HERE, 'results_eta.csv')).read())

    hdr = ('| L | β | η̂ (row sum) | η max | **witness lb** (max link) | cert/est '
           '| mean TV/pair | mean ‖ΔS‖_* |\n'
           '|---|---|---|---|---|---|---|---|')
    lines = [hdr]
    for r in rows:
        lines.append(f"| {r['L']} | {r['beta']:.1f} | {r['eta']:.2f}({r['eta_err']:.2f}) | "
                     f"{r['eta_max']:.2f} | **{r['lb_max']:.2f}** | {r['cert']} | "
                     f"{r['tv_mean']:.3f} | {r['ds_mean']:.2f} |")
    txt = '\n'.join(lines)
    open(os.path.join(HERE, 'na3t_table.txt'), 'w').write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
