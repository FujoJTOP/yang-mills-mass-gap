"""Negative-mode statistics of the SU(3) Wilson Hessian at the crossover.

The MC samples the Gibbs measure, so ensemble averages ARE measure-weighted.
Reports, per link (link_hessian_min = lowest eigenvalue of the 8x8 per-link
Hessian): mean/std/min/max, fraction with min<0, E[(-x)_+] (negative part),
E[(x)_+] (positive part), their ratio, and the per-config mean.
"""
import sys

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import numpy as np  # noqa: E402
import su3 as S  # noqa: E402


def run(beta, L, sweeps, therm, seed=11):
    rng = np.random.default_rng(seed)
    lat = S.Lattice(L, beta, rng, cold=True)
    for _ in range(therm):
        lat.sweep()
    mins = []
    for _ in range(sweeps):
        lat.sweep()
        mins.append(S.link_hessian_min(lat.U, beta).ravel().astype(float))
    arr = np.stack(mins)                    # (sweeps, nlinks)
    m = arr.ravel()
    neg = np.clip(-m, 0.0, None)
    pos = np.clip(m, 0.0, None)
    cfg_mean = arr.mean(axis=1)
    print(f"beta={beta} L={L}  nlinks={arr.shape[1]} samples={m.size}")
    print(f"  mean={m.mean():+.4f} std={m.std():.4f} min={m.min():+.4f} max={m.max():+.4f}")
    print(f"  frac(min<0)={np.mean(m < 0):.4f}  E[neg]={neg.mean():.4f}  "
          f"E[pos]={pos.mean():.4f}  E[neg]/E[pos]={neg.mean()/max(pos.mean(),1e-12):.3f}")
    print(f"  q01={np.quantile(m,0.01):+.4f} q05={np.quantile(m,0.05):+.4f} "
          f"q50={np.quantile(m,0.5):+.4f}")
    print(f"  per-config mean: {cfg_mean.mean():+.4f}, "
          f"frac(config mean<0)={np.mean(cfg_mean < 0):.4f}")
    return dict(beta=beta, frac_neg=float(np.mean(m < 0)),
                Eneg=float(neg.mean()), Epos=float(pos.mean()),
                min=float(m.min()), mean=float(m.mean()))


if __name__ == '__main__':
    for b in (2.5, 2.6, 3.0):
        run(b, 6, 80, 300)
