"""Focused final analysis on a few representative modes, with data saved and a
conditional-window scan that verifies  kappa_hist  ==  E[H|a~0] - Var(G|a~0).

For each mode:
  kappa_hist, SE (bootstrap), lam4 (quartic), 1/Var(a);
  E[H_vv], Var(G_v)  (bare, log-free);
  width-scan:  Ev(w) = mean(H | |a|<w*sd) - var(G | |a|<w*sd)  for several w,
               -> should approach kappa_hist as w->0 (confirms identity+binning).
"""
import sys
import time
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import w4_blocked as W    # noqa: E402
import gpp_core as C      # noqa: E402
import gpp_core3 as G3    # noqa: E402
import w7_run3 as R3      # noqa: E402


def width_scan(a, H, G, ws=(0.1, 0.15, 0.2, 0.3, 0.5)):
    sd = a.std()
    out = []
    for w in ws:
        m = np.abs(a) < w * sd
        if m.sum() < 30:
            out.append((w, np.nan, np.nan, int(m.sum())))
            continue
        out.append((w, float(H[m].mean()), float(G[m].var()), int(m.sum())))
    return out


def analyse_mode(a, H, G, rng):
    kap, se, lam, mb, inv = C.veff_curvature(a, rng=rng, nboot=500)
    scan = width_scan(a, H, G)
    return dict(kap=kap, se=se, lam=lam, inv=inv, EH=H.mean(), VG=G.var(),
                sd=a.std(), scan=scan, Ev_uncond=H.mean() - G.var())


SU2_MODES = {
    'ax1_p2_c2_cos': (1, 2, 2, 'cos'),
    'ax1_p3_c2_cos': (1, 3, 2, 'cos'),
    'ax0_p2_c2_cos': (0, 2, 2, 'cos'),
    'ax1_p0_c2_cos': (1, 0, 2, 'cos'),
    'ax1_p2_c0_cos': (1, 2, 0, 'cos'),
    'ax2_p3_c2_sin': (2, 3, 2, 'sin'),
}


def run_su2(L=4, beta=2.5, N=4000, seed=303):
    cfgs = C.sample_configs(L, beta, 400, N, 3, seed=seed)
    names = list(SU2_MODES)
    ml = [C.plane_mode(L, *SU2_MODES[n][:3], n=1, kind=SU2_MODES[n][3]) for n in names]
    acc = {n: ([], [], []) for n in names}
    t0 = time.time()
    for U in cfgs:
        a, H, G = C.mode_obs(U, beta, ml, h=1e-3)
        for n, av, hv, gv in zip(names, a, H, G):
            acc[n][0].append(av); acc[n][1].append(hv); acc[n][2].append(gv)
    print(f"SU2 done in {time.time()-t0:.0f}s", flush=True)
    rng = np.random.default_rng(9001)
    np.savez('su2_data.npz', **{f'{n}_{k}': np.array(acc[n][i])
              for n in names for i, k in enumerate(('a', 'H', 'G'))})
    for n in names:
        a, H, G = map(np.array, acc[n])
        r = analyse_mode(a, H, G, rng)
        print(f"\n[{n}]  Gpp={r['kap']:+.4f}+-{r['se']:.4f}  1/Var={r['inv']:+.4f} "
              f"lam4={r['lam']:+.4f}  E[H]={r['EH']:.4f} VarG={r['VG']:.4f} sd={r['sd']:.3f}")
        print(f"   width-scan Ev(w)=E[H|w]-Var(G|w): " +
              "  ".join(f"w={w:.2f}:{ev:+.3f}({cnt})" if not np.isnan(ev)
                        else f"w={w}:nan" for (w, eh, vg, cnt) in r['scan']
                        for ev in [eh - vg]))


def run_su3(L=4, beta=5.7, N=1500, seed=404):
    cfgs = G3.sample_cfgs(L, beta, 200, N, 3, seed=seed)
    spec = {'ax1_p2_c3_cos': (1, 2, 3), 'ax1_p3_c3_cos': (1, 3, 3),
            'ax0_p2_c3_cos': (0, 2, 3), 'ax1_p0_c3_cos': (1, 0, 3),
            'ax1_p2_c0_cos': (1, 2, 0), 'ax1_p2_c7_cos': (1, 2, 7)}
    names = list(spec)
    ml = [G3.plane_mode(L, *spec[n], n=1, kind='cos') for n in names]
    acc = {n: ([], [], []) for n in names}
    ang = []
    t0 = time.time()
    for U in cfgs:
        ang.append(R3.max_eigangle(U))
        a, H, G = G3.mode_obs(U, beta, ml, h=1e-3)
        for n, av, hv, gv in zip(names, a, H, G):
            acc[n][0].append(av); acc[n][1].append(hv); acc[n][2].append(gv)
    ang = np.array(ang); safe = ang < (np.pi - 0.15)
    print(f"SU3 done {time.time()-t0:.0f}s; branch-safe {safe.sum()}/{N}", flush=True)
    rng = np.random.default_rng(9002)
    np.savez('su3_data.npz', **{f'{n}_{k}': np.array(acc[n][i])
              for n in names for i, k in enumerate(('a', 'H', 'G'))})
    for n in names:
        a, H, G = map(np.array, acc[n])
        r = analyse_mode(a, H, G, rng)
        rs = analyse_mode(a[safe], H[safe], G[safe], rng) if safe.sum() > 200 else None
        print(f"\n[{n}]  Gpp={r['kap']:+.4f}+-{r['se']:.4f}  1/Var={r['inv']:+.4f} "
              f"E[H]={r['EH']:.4f} VarG={r['VG']:.4f} sd={r['sd']:.3f}"
              + (f"  | safe: {rs['kap']:+.4f}+-{rs['se']:.4f}" if rs else ""))


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--group', default='su2')
    ap.add_argument('--N', type=int, default=0)
    ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    if a.group == 'su2':
        run_su2(N=a.N or 4000, seed=a.seed or 303)
    else:
        run_su3(N=a.N or 1500, seed=a.seed or 404)
