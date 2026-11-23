"""First measurement for the local-Gaussianity entry (ENTRY.md section 4).

At the crossover (SU(3), beta = 2.6 / 3.0 / 3.3) measure:
  * per-link Hessian lowest eigenvalue distribution  -> how deep / how often the
    Savvidy negative modes are present in the *microscopic* action;
  * blocked plaquette excess kurtosis R4(b)           -> Gaussianity of the
    blocked field (b = 1, 2);
  * force (dS_W) excess kurtosis                      -> Gaussianity of the
    bare force.

Honest scope: the Hessian min is the MICROSCOPIC Hessian (not the blocked
action's).  R4 is the blocked plaquette (a gauge-invariant slice of the blocked
field).  Neither is the full blocked Hessian; together they are a first,
runnable probe of the "effectively SLC / Gaussian" hypothesis.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'lattice-cpb'))
import su3 as S  # noqa: E402


def block_mean(P, b):
    L = P.shape[0]
    n = L // b
    return P.reshape(n, b, n, b, n, b, n, b).mean(axis=(1, 3, 5, 7))


def excess_kurtosis(x):
    x = np.asarray(x, dtype=float).ravel()
    m = x.mean()
    v = np.mean((x - m) ** 2)
    if v <= 0:
        return 0.0
    return float(np.mean((x - m) ** 4) / v ** 2 - 3.0)


def run(beta, L, sweeps, therm, seed):
    rng = np.random.default_rng(seed)
    lat = S.Lattice(L, beta, rng, cold=True)
    for _ in range(therm):
        lat.sweep()

    hmin_all, plaq, Pb1, Pb2, force_k4 = [], [], [], [], []
    for _ in range(sweeps):
        lat.sweep()
        U = lat.U
        hmin = S.link_hessian_min(U, beta)
        hmin_all.append(hmin.ravel())
        tr = S.plaquette_traces(U)
        P = np.mean(np.stack(list(tr.values())), axis=0)
        plaq.append(P.mean())
        Pb1.append(P.copy())
        Pb2.append(block_mean(P, 2).copy())
        G = S.action_gradient(U, beta)          # (4,L,L,L,L,8)
        force_k4.append(G.reshape(-1, G.shape[-1]))

    hmin_all = np.concatenate(hmin_all)          # flat over links*configs
    fk = np.concatenate(force_k4, axis=0)        # (N,8)
    # excess kurtosis of each force component
    fk4 = [excess_kurtosis(fk[:, a]) for a in range(fk.shape[1])]

    # blocked plaquette: per-config excess kurtosis over sites, averaged over configs
    def r4_mean(lst):
        vals = []
        for arr in lst:
            vals.append(excess_kurtosis(arr))
        return float(np.mean(vals))

    return {
        'beta': beta,
        'plaq': float(np.mean(plaq)),
        'hmin_mean': float(hmin_all.mean()),
        'hmin_std': float(hmin_all.std()),
        'hmin_min': float(hmin_all.min()),
        'hmin_frac_neg': float(np.mean(hmin_all < 0.0)),
        'hmin_p01': float(np.quantile(hmin_all, 0.01)),
        'hmin_p99': float(np.quantile(hmin_all, 0.99)),
        'R4_b1': r4_mean(Pb1),
        'R4_b2': r4_mean(Pb2),
        'force_k4_mean': float(np.mean(fk4)),
        'force_k4_absmax': float(np.max(np.abs(fk4))),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--L', type=int, default=4)
    ap.add_argument('--sweeps', type=int, default=400)
    ap.add_argument('--therm', type=int, default=800)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--betas', type=str, default='2.6,3.0,3.3')
    args = ap.parse_args()

    print(f"L={args.L} sweeps={args.sweeps} therm={args.therm}")
    hdr = ("beta   <P>      hmin_mean  hmin_std  hmin_min  frac<0  p01     "
           "R4_b1     R4_b2     force_k4")
    print(hdr)
    for b in [float(x) for x in args.betas.split(',')]:
        r = run(b, args.L, args.sweeps, args.therm, args.seed)
        print(f"{r['beta']:5.2f}  {r['plaq']:.4f}   {r['hmin_mean']:+.4f}   "
              f"{r['hmin_std']:.4f}   {r['hmin_min']:+.4f}   "
              f"{r['hmin_frac_neg']:.3f}   {r['hmin_p01']:+.4f}   "
              f"{r['R4_b1']:+.4f}   {r['R4_b2']:+.4f}   {r['force_k4_mean']:+.3f}")
        print(f"      force k4 absmax = {r['force_k4_absmax']:.3f}", flush=True)


if __name__ == '__main__':
    main()
