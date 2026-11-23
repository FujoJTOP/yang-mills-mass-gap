"""
RERUN driver: measure G''(0) and its decomposition for SU(2)/SU(3), L=4/6.

Per configuration and mode v (flat exponential-chart coordinates):
  a      = v . x                     (x = log-algebra coords, principal branch)
  H_vv   = flat  d^2/ds^2 S(X+sV)    (action Hessian)
  G_v    = flat  d/ds    S(X+sV)    (action force)
  Phi'   = flat  d/ds    (-logJ)(X+sV)
  Phi''  = flat  d^2/ds^2(-logJ)(X+sV)
Then:  V_eff''(0) = E[H+Phi''|a=0] - Var(G+Phi'|a=0)
                  = (E[H]-Var(G)) + W ,   W = E[Phi'']-Var(Phi')-2Cov(G,Phi')
Gold standard: kappa_hist = curvature of -log P(a).
"""
import sys
import time
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import w4_blocked as W        # noqa: E402
import gpp_core as C          # noqa: E402
import gpp_core3 as G3        # noqa: E402
import su3 as S3              # noqa: E402
import rerun_core as RC       # noqa: E402


# ======================================================================
# sampling
# ======================================================================
def sample_su2(L, beta, ntherm, ncfg, gap, seed):
    rng = np.random.default_rng(seed)
    return W.sample2(L, beta, rng, ntherm, ncfg, gap)


def sample_su3(L, beta, ntherm, ncfg, gap, seed, cold=False):
    rng = np.random.default_rng(seed)
    lat = S3.Lattice(L, beta, rng, cold=cold)
    for _ in range(ntherm):
        lat.sweep()
    out = []
    for _ in range(ncfg):
        for _ in range(gap):
            lat.sweep()
        out.append(lat.U.copy())
    return out


# ======================================================================
# per-config measurement for a set of modes
# ======================================================================
def measure(U, beta, modes, group, h=1e-3):
    """returns dict of arrays [a,H,G,P1,P2] per mode, plus branch diagnostics."""
    if group == 'su2':
        x = C.su2_log_x(U)
        Xfrom = RC.X_from_x_su2
        action = W.action2
    else:
        _, x = RC.su3_log_herm(U)          # branch-safe traceless log
        Xfrom = RC.X_from_x_su3
        action = G3.S3.wilson_action
    out = {k: [] for k in ('a', 'H', 'G', 'P1', 'P2')}
    if group == 'su3':
        w = np.linalg.eigvals(U)
        s = np.round(np.angle(w).sum(axis=-1) / (2 * np.pi))
        out['ntwist'] = int(np.count_nonzero(s))
        out['nlink'] = int(s.size)
    for v in modes:
        a = float(np.sum(v * x))
        G, H = RC.flat_derivs(x, v, beta, Xfrom, action, h)
        P1, P2 = RC.phi_derivs_analytic(x, v, Xfrom)
        out['a'].append(a); out['H'].append(H); out['G'].append(G)
        out['P1'].append(P1); out['P2'].append(P2)
    return {k: np.array(val) for k, val in out.items()}


# ======================================================================
# analysis: gold-standard kappa + decomposition, conditioned on a~0
# ======================================================================
def analyse_mode(a, H, G, P1, P2, rng, bw=0.30):
    kap, se, lam4, inv = RC.kappa_hist(a, rng=rng)
    kkde, skde = RC.kappa_kde(a, rng=rng)
    sd = a.std()
    # kernel-conditioned moments at a=0 (the exact identity needs a=0)
    k = RC.kernel_cond(a, {'H': H, 'G': G, 'P1': P1, 'P2': P2}, bw=bw)
    EH = k['E_H']; EP2 = k['E_P2']
    VG = k['V_own']['G']; VP1 = k['V_own']['P1']; cov = k['C_G_P1']
    bracket = EH - VG                 # action part at a=0
    W = EP2 - VP1 - 2.0 * cov         # Jacobian/geometric part at a=0
    tot = (EH + EP2) - (VG + VP1 + 2.0 * cov)   # = bracket + W  (identity RHS)
    cond = dict(EH=EH, VG=VG, EP2=EP2, VP1=VP1, cov=cov,
                bracket=bracket, W=W, tot=tot, n_eff=k['n_eff'])
    # unconditioned (for comparison / to show why conditioning matters)
    uH, uG, uP1, uP2 = H.mean(), G.var(), float(np.var(P1)), P2.mean()
    ucov = float(np.mean((G - G.mean()) * (P1 - P1.mean())))
    unc = dict(EH=uH, VG=uG, EP2=uP2, VP1=uP1, cov=ucov,
               bracket=uH - uG, W=uP2 - uP1 - 2 * ucov,
               tot=(uH + uP2) - (uG + uP1 + 2 * ucov))
    return dict(kap=kap, se=se, lam4=lam4, inv=inv, sd=sd,
                kkde=kkde, skde=skde, cond=cond, unc=unc)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--group', default='su2')
    ap.add_argument('--L', type=int, default=4)
    ap.add_argument('--beta', type=float, default=None)
    ap.add_argument('--N', type=int, default=600)
    ap.add_argument('--therm', type=int, default=400)
    ap.add_argument('--gap', type=int, default=3)
    ap.add_argument('--seed', type=int, default=101)
    args = ap.parse_args()
    if args.group == 'su2':
        beta = 2.5 if args.beta is None else args.beta
        cfgs = sample_su2(args.L, beta, args.therm, args.N, args.gap, args.seed)
        # 6 representative transverse modes
        spec = [(1, 2, 2, 'cos'), (1, 3, 2, 'cos'), (0, 2, 2, 'cos'),
                (2, 3, 2, 'cos'), (3, 1, 2, 'cos'), (1, 2, 0, 'cos')]
        modes = [C.plane_mode(args.L, ax, pl, c, n=1, kind=k) for (ax, pl, c, k) in spec]
        names = [f'ax{ax}_p{pl}_c{c}_{k}' for (ax, pl, c, k) in spec]
    else:
        beta = 5.7 if args.beta is None else args.beta
        cfgs = sample_su3(args.L, beta, args.therm, args.N, args.gap, args.seed)
        spec = [(1, 2, 3), (1, 3, 3), (0, 2, 3), (2, 3, 7), (3, 1, 3), (1, 2, 7)]
        modes = [G3.plane_mode(args.L, ax, pl, c, n=1, kind='cos') for (ax, pl, c) in spec]
        names = [f'ax{ax}_p{pl}_c{c}' for (ax, pl, c) in spec]
    print(f"# {args.group} L={args.L} beta={beta} N={args.N} modes={len(modes)}", flush=True)
    print(f"  sampled {len(cfgs)} cfgs", flush=True)
    t0 = time.time()
    acc = {k: [] for k in ('a', 'H', 'G', 'P1', 'P2')}
    tw = []
    for i, U in enumerate(cfgs):
        r = measure(U, beta, modes, args.group)
        for k in acc:
            acc[k].append(r[k])
        if args.group == 'su3':
            tw.append(r['ntwist'])
        if (i + 1) % 100 == 0:
            print(f"    {i+1}/{args.N} ({time.time()-t0:.0f}s)", flush=True)
    if args.group == 'su3':
        tw = np.array(tw)
        print(f"  center-twist: {tw.sum()}/{len(tw)*r['nlink']} links "
              f"({100*tw.sum()/(len(tw)*r['nlink']):.2f}%); configs with >=1 twist: "
              f"{np.count_nonzero(tw)}/{len(tw)}; corrected log is exact (expi(X)=U) for ALL",
              flush=True)
    # stack: (Ncfg, Nmode)
    for k in acc:
        acc[k] = np.array(acc[k])
    np.savez(f'rerun_{args.group}_L{args.L}.npz', **acc)
    rng = np.random.default_rng(7000)
    print(f"\n  {'mode':22s} {'Gpp(kde)':>9s} {'+-SE':>6s} {'kap_hist':>8s} {'1/Var':>6s} "
          f"{'E[H]':>6s} {'VarG':>6s} {'EP2':>6s} {'br':>7s} {'W':>7s} {'sum_c':>7s}", flush=True)
    for j, nm in enumerate(names):
        a = acc['a'][:, j]; H = acc['H'][:, j]; G = acc['G'][:, j]
        P1 = acc['P1'][:, j]; P2 = acc['P2'][:, j]
        r = analyse_mode(a, H, G, P1, P2, rng)
        c = r['cond']
        print(f"  {nm:22s} {r['kkde']:+9.4f} {r['skde']:6.4f} {r['kap']:+8.4f} {r['inv']:6.3f} "
              f"{c['EH']:6.3f} {c['VG']:6.3f} {c['EP2']:6.3f} {c['bracket']:+7.3f} "
              f"{c['W']:+7.3f} {c['tot']:+7.3f}", flush=True)
