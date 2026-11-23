"""Weak-coupling control for kappa_eff: is kappa_eff > 0 where the theory is
known to be fine?  (beta = 8,10,20,40, L=6)"""
import numpy as np

import su3 as S
import kappa as K

L = 6
recs = []
for beta in (8.0, 10.0, 20.0, 40.0):
    rng = np.random.default_rng(int(1000 * beta))
    lat = S.Lattice(L, beta, rng, cold=True)
    for _ in range(600):
        lat.sweep()
    U = lat.U.copy()
    hl = S.link_hessian_min(U, beta)
    res = K.lowest_eigenvalue(U, beta, m=160, seed=1)
    rec = dict(L=L, beta=beta, snap=0, kappa=res['kappa'],
               lanczos_res=res['lanczos_res'], cg_res=res['cg_res_max'],
               ritz=res['ritz'], hl_min=float(hl.min()),
               hl_neg_frac=float(np.mean(hl < 0)), hl_mean=float(hl.mean()),
               plaq=S.mean_plaquette(U))
    recs.append(rec)
    print(f"beta={beta}: P={rec['plaq']:.5f} kappa={rec['kappa']:+.5f} "
          f"res={rec['lanczos_res']:.1e} linkH mean={rec['hl_mean']:+.3f} "
          f"fracneg={rec['hl_neg_frac']:.3f}", flush=True)

out = {}
for key in ('L', 'beta', 'snap', 'kappa', 'lanczos_res', 'cg_res',
            'hl_min', 'hl_neg_frac', 'hl_mean'):
    out[key] = np.array([r[key] for r in recs])
out['ritz'] = np.array([r['ritz'] for r in recs])
out['plaq'] = np.array([r['plaq'] for r in recs])
np.savez('kappa_results_control.npz', **out)
print('saved kappa_results_control.npz')
