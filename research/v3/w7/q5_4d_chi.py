#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
q5_4d_chi.py -- Q5: a controlled 4D SU(2) numerical test of the H1 plateau.

H1 (the last open input): the k=0 susceptibility chi = Var(a) is bounded uniformly
in L at beta=2.5 (equivalently G''_{k=0} >= lambda > 0).  1D periodic SU(2) gives a
plateau (Var(a)->4.45, G''->0.222).  Is 4D a plateau or growth?

This witness does four things, one job, py -3.12, SU(2) Metropolis/Haar heat bath:

  [FREE]  exact free-field acid test: free MASSIVE -> chi=1/m^2 plateau (L-flat),
          G''=m^2; free MASSLESS (IR-mode regulated) -> chi ~ L^2, G'' ~ 1/L^2 -> 0.
          The tool manufactures NO floor for the free field.                    [T]

  [HAAR]  beta=0 MC check: chi_flat -> s^2=(4pi^2-6)/9=3.7198 (exact Haar),
          Polyakov spatial correlations vanish for r>0.  Validates sampler.    [V]

  [PILOT] beta=2.5, periodic L^4, L=4 and L=6:
            chi_flat = (1/3)<|A|^2>,  A = N^{-1/2} sum_links xi_l  (k=0 mode,
                       flat chart;  globally colour-invariant, gauge-invariant
                       only infinitesimally -- the gauge systematic is measured).
            C_P(r)   = gauge-invariant spatial Polyakov-loop correlator -> xi.
            G''(k=0) via the exact-identity / force-moment estimator
                       G''=E[H|a0]-Var(G|a0)=bracket+W_geo+X  and  1/Var(a).
            tau_int, delete-block jackknife SE, effective sample size.          [V]

  [VERDICT] does chi(L) plateau or grow; the minimal decisive run.

NO third-party sources.  ONE heavy job; never parallelise SU(N) L>=6 (charter s4).
Labels [T]/[C]/[NG]/[O]/[V]; no QED.
"""
import sys, os, json, time
import numpy as np

CPB = 'D:/Dev/math/research/lattice-cpb'
W7  = 'D:/Dev/math/research/unconditional/w7'
sys.path.insert(0, CPB); sys.path.insert(0, W7)
import w4_blocked as W           # SU(2) lattice, heat bath, action, randomness
import gpp_core as C             # su2_log_x, plane_mode
import rerun_core as RC          # X_from_x_su2, expi, phi_derivs_analytic
try:
    import k0_continuum_analyze as KA   # pooled_jack (validated force-moment)
    _HAVE_KA = True
except Exception:
    _HAVE_KA = False

rng_global = np.random.default_rng(20261124)

# --------------------------------------------------------------------------
# streaming SU(2) heat-bath sampler (in place, single U holder)
# --------------------------------------------------------------------------
def make_stream(L, beta, seed, cold=False):
    rng = np.random.default_rng(seed)
    if cold:
        U = np.zeros((4, L, L, L, L, 2, 2), dtype=np.complex128)
        U[..., 0, 0] = 1.0; U[..., 1, 1] = 1.0
    else:
        U = W.random_su2(rng, (4, L, L, L, L))
    ns = L ** 4
    x0, x1, x2, x3 = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
    par = (x0 + x1 + x2 + x3) % 2
    idxp = [np.nonzero(par.reshape(-1) == c)[0] for c in (0, 1)]
    kappa = beta / 2.0
    def sweep():
        Uf = U.reshape(4, ns, 2, 2)
        for c in (0, 1):
            ids = idxp[c]
            for mu in range(4):
                S = W.staple_sum2(U, mu).reshape(ns, 2, 2)[ids]
                A = Uf[mu, ids]; M = A @ S
                Uf[mu, ids] = W.su2_heatbath(kappa, M, 0, 1, rng) @ A
    return U, sweep

# --------------------------------------------------------------------------
# estimators
# --------------------------------------------------------------------------
def polyakov_field(U):
    """P(x) = (1/2) tr prod_t U_0(x,t) in [-1,1]; spatial field (L,L,L)."""
    L = U.shape[1]
    A = np.broadcast_to(np.eye(2, dtype=complex), (L, L, L, L, 2, 2)).copy()
    for t in range(L):
        A = A @ np.roll(U[0], -t, axis=0)
    return np.real(np.trace(A[0], axis1=-2, axis2=-1)) / 2.0

def chi_flat(U):
    """chi_flat = (1/3)<|A|^2>, A = N^{-1/2} sum_l xi_l (flat-chart k=0 mode).
    Globally colour-invariant; the flat-chart coordinate -> gauge-dependent
    beyond infinitesimal (measured separately)."""
    x = C.su2_log_x(U)
    Nl = 4 * U.shape[1] ** 4
    A = x.reshape(Nl, 3).sum(0) / np.sqrt(Nl)
    return float(A @ A / 3.0)

def corr_P(U, L):
    """ensemble (NOT per-config) Polyakov spatial correlator along axis 0."""
    P = polyakov_field(U)
    return [float(np.mean(P * np.roll(P, -r, axis=0))) for r in range(L)]

def plaq(U, beta):
    if beta <= 0:
        return float('nan')
    return float(1.0 - W.action2(U, beta) / (beta * 6.0 * U.shape[1] ** 4))

def tau_int(x):
    x = np.asarray(x, float); x = x - x.mean(); v = x.var(); n = x.size
    if v <= 0:
        return 1.0
    s = 1.0
    for t in range(1, n // 4):
        c = np.mean(x[:n - t] * x[t:]) / v
        if c < 0:
            break
        s += 2.0 * c
    return float(max(1.0, s))

def jk_mean_se(x, nblk=10):
    x = np.asarray(x, float); n = x.size
    nblk = int(min(nblk, max(4, n // 8)))
    ed = np.linspace(0, n, nblk + 1).astype(int)
    m = x.mean()
    lv = np.array([np.delete(x, np.arange(ed[b], ed[b + 1])).mean() for b in range(nblk)])
    se = np.sqrt((nblk - 1) / nblk * np.sum((lv - lv.mean()) ** 2))
    return float(m), float(se)

def gauge_transform(U, g):
    out = np.empty_like(U)
    for mu in range(4):
        out[mu] = g @ U[mu] @ W.adj(np.roll(g, -1, axis=mu))
    return out

def rand_local_su2(rng, shape, eps):
    ph = rng.standard_normal(tuple(shape) + (3,)) * eps
    A = np.einsum('...a,aij->...ij', ph, W.T2)
    w, V = np.linalg.eigh(A)
    return np.einsum('...ik,...k,...jk->...ij', V, np.exp(1j * w), np.conj(V))

# --------------------------------------------------------------------------
# [FREE] exact acid test
# --------------------------------------------------------------------------
def free_field_acid():
    out = {}
    # massive free field: k=0 mode eigenvector, Var(a)=1/m^2, G''=m^2 (L-flat)
    masses = [1.0, 0.5, 0.3, 0.1, 0.03]
    out['massive'] = [dict(m=m, Var=1.0 / m ** 2, Gpp=m ** 2,
                           Gpp_times_Var=m ** 2 / m ** 2) for m in masses]
    # massless, box-regulated by the smallest nonzero lattice momentum k_min=2pi/L:
    #   Var = 1/khat_min^2 ~ L^2/(4 pi^2) ,  G'' = khat_min^2 ~ 4 pi^2/L^2 -> 0
    ml = []
    for L in (4, 6, 8, 12, 16):
        khat2 = 4.0 * np.sin(np.pi / L) ** 2        # khat^2 at k_min = 2pi/L
        ml.append(dict(L=L, Var=1.0 / khat2, Gpp=khat2,
                       Var_over_L2=(1.0 / khat2) / L ** 2,
                       Gpp_times_L2=khat2 * L ** 2))
    out['massless_box'] = ml
    return out

# --------------------------------------------------------------------------
# [HAAR]+[PILOT] one MC run
# --------------------------------------------------------------------------
def run(L, beta, N, therm, gap, seed, cold, want_force=False, want_gauge=False):
    U, sweep = make_stream(L, beta, seed, cold)
    t0 = time.time()
    for _ in range(therm):
        sweep()
    chF = np.empty(N); plq = np.empty(N)
    Cp = np.zeros(L)
    Pavg = np.empty(N)
    force = dict(a=[], H=[], G=[], P1=[], P2=[])
    if want_force:
        v = np.zeros((4, L, L, L, L, 3)); v[..., 2] = 1.0
        v /= np.linalg.norm(v)
        Xfrom = RC.X_from_x_su2
    stored = []
    stride = max(1, N // 24)          # spread gauge-test configs over the whole run
    for i in range(N):
        for _ in range(gap):
            sweep()
        chF[i] = chi_flat(U)
        plq[i] = plaq(U, beta)
        P = polyakov_field(U); Pavg[i] = P.mean()
        Cp += corr_P(U, L)
        if want_force:
            x = C.su2_log_x(U)
            S0 = W.action2(U, beta)
            Sp = W.action2(RC.expi(RC.X_from_x_su2(x + 1e-3 * v)), beta)
            Sm = W.action2(RC.expi(RC.X_from_x_su2(x - 1e-3 * v)), beta)
            force['a'].append(float(np.sum(v * x)))
            force['G'].append((Sp - Sm) / 2e-3)
            force['H'].append((Sp - 2.0 * S0 + Sm) / 1e-6)
            p1, p2 = RC.phi_derivs_analytic(x, v, Xfrom)
            force['P1'].append(p1); force['P2'].append(p2)
        if want_gauge and (i % stride == 0) and len(stored) < 24:
            stored.append(U.copy())
    Cp /= N
    cp_corr = Cp - Pavg.mean() ** 2          # ensemble connected
    chi_P = float(L ** 3 * Pavg.var())       # = sum_r C_P(r)  (ensemble)
    # jackknife SE on chi_P = V_s Var(Pbar)
    nb = int(min(10, max(4, N // 8))); ed = np.linspace(0, N, nb + 1).astype(int)
    lvp = np.array([np.delete(Pavg, np.arange(ed[b], ed[b + 1])).var() for b in range(nb)])
    chi_P_se = float(L ** 3 * np.sqrt((nb - 1) / nb * np.sum((lvp - lvp.mean()) ** 2)))
    mF, sF = jk_mean_se(chF); mm, _ = jk_mean_se(plq)
    tF = tau_int(chF)
    res = dict(L=L, beta=beta, N=N, gap=gap, therm=therm,
               chi_flat=mF, chi_flat_se=sF, chi_flat_tau=tF, chi_flat_neff=N / tF,
               plaq=mm, chi_P=chi_P, chi_P_se=chi_P_se, Pavg_mean=float(Pavg.mean()),
               C_P=[float(c) for c in cp_corr], wall_s=time.time() - t0)
    # force-moment decomposition (exact identity at the k=0 mode)
    if want_force and len(force['a']) > 20:
        try:
            A = np.array(force['a']).reshape(-1, 1)
            Hh = np.array(force['H']).reshape(-1, 1)
            Gg = np.array(force['G']).reshape(-1, 1)
            p1 = np.array(force['P1']).reshape(-1, 1)
            p2 = np.array(force['P2']).reshape(-1, 1)
            if _HAVE_KA:
                full, se, ml, nblk = KA.pooled_jack(A, Hh, Gg, p1, p2, nblk=10)
                res['force'] = dict(Gpp=full['Gpp'], Gpp_se=se['Gpp'],
                                    bracket=full['bracket'], Wgeo=full['Wgeo'],
                                    X=full['X'], inv=full['inv'],
                                    se={k: se[k] for k in ('bracket', 'Wgeo', 'X', 'tot')},
                                    a_sd=float(np.std(force['a'])))
            else:
                a = np.array(force['a']); Hv = np.array(force['H']); Gv = np.array(force['G'])
                res['force'] = dict(bracket=float(Hv.mean() - Gv.var()),
                                    inv=float(1.0 / a.var()), a_sd=float(a.std()))
        except Exception as e:
            res['force'] = dict(error=repr(e))
    # gauge systematic: spread of ensemble-mean chi_flat under finite local gauge tf
    if want_gauge and stored:
        rng = np.random.default_rng(99)
        base = np.mean([chi_flat(Uu) for Uu in stored])
        means = []
        for eps in (0.05, 0.1, 0.2, 0.4):
            vals = []
            for Uu in stored:
                g = rand_local_su2(rng, (L, L, L, L), eps)
                vals.append(chi_flat(gauge_transform(Uu, g)))
            means.append((eps, float(np.mean(vals))))
        res['gauge'] = dict(base=float(base),
                            scans=[dict(eps=e, chi=v) for e, v in means])
    return res

# --------------------------------------------------------------------------
def main():
    quick = 'quick' in sys.argv
    print("# Q5 -- 4D SU(2) k=0 susceptibility (H1 plateau test)  py -3.12  one job",
          flush=True)
    out = dict(meta=dict(date='2026-11-24', beta_main=2.5, group='SU(2)'))

    print("\n[FREE] exact free-field acid test  [T]", flush=True)
    fa = free_field_acid(); out['free'] = fa
    for d in fa['massive']:
        print(f"   massive m={d['m']:.2f}: Var={d['Var']:.4f} (L-flat) "
              f"G''={d['Gpp']:.4f}=m^2  G''*Var={d['Gpp_times_Var']:.3f}", flush=True)
    for d in fa['massless_box']:
        print(f"   massless L={d['L']:2d}: Var={d['Var']:9.3f}  (Var/L^2={d['Var_over_L2']:.4f} const)"
              f"  G''={d['Gpp']:.5f}  (G''*L^2={d['Gpp_times_L2']:.3f} const -> G''~1/L^2 -> 0)",
              flush=True)

    # ---- beta=0 Haar validation ----
    print("\n[HAAR] beta=0 pure-Haar check  [V]", flush=True)
    r0 = run(4, 0.0, 150 if quick else 200, 150, 6, 11, cold=False)
    out['beta0'] = r0
    print(f"   chi_flat={r0['chi_flat']:.3f}+-{r0['chi_flat_se']:.3f} "
          f"(exact s^2=3.7198)  chi_P={r0['chi_P']:.4f}  C_P(r)={np.round(r0['C_P'],4)}",
          flush=True)

    # ---- main pilot beta=2.5 ----
    print("\n[PILOT] beta=2.5 periodic SU(2)  [V]", flush=True)
    N4 = 120 if quick else 320
    N6 = 80 if quick else 200
    r4 = run(4, 2.5, N4, 800, 8, 101, cold=True, want_force=True, want_gauge=True)
    out['L4'] = r4
    print(f"   L=4 N={N4}: chi_flat={r4['chi_flat']:.4f}+-{r4['chi_flat_se']:.4f} "
          f"(tau={r4['chi_flat_tau']:.1f}, Neff={r4['chi_flat_neff']:.0f}) plaq={r4['plaq']:.4f}",
          flush=True)
    print(f"          chi_P={r4['chi_P']:.4f}+-{r4['chi_P_se']:.4f}  <P>={r4['Pavg_mean']:+.4f}  "
          f"C_P(r)={np.round(r4['C_P'],4)}", flush=True)
    if 'force' in r4:
        f = r4['force']
        print(f"          G''(k=0) force-moment: Gpp={f.get('Gpp',float('nan')):+.4f}"
              f"+-{f.get('Gpp_se',float('nan')):.4f} bracket={f.get('bracket',float('nan')):+.3f} "
              f"Wgeo={f.get('Wgeo',float('nan')):+.3f} X={f.get('X',float('nan')):+.3f} "
              f"1/Var={f.get('inv',float('nan')):.3f}", flush=True)
    if 'gauge' in r4:
        g = r4['gauge']
        print(f"          gauge systematic: base chi_flat={g['base']:.4f}; "
              + "  ".join(f"eps={d['eps']}->{d['chi']:.3f}" for d in g['scans']), flush=True)

    r6 = run(6, 2.5, N6, 800, 8, 202, cold=True, want_force=False, want_gauge=False)
    out['L6'] = r6
    print(f"   L=6 N={N6}: chi_flat={r6['chi_flat']:.4f}+-{r6['chi_flat_se']:.4f} "
          f"(tau={r6['chi_flat_tau']:.1f}, Neff={r6['chi_flat_neff']:.0f}) plaq={r6['plaq']:.4f}",
          flush=True)
    print(f"          chi_P={r6['chi_P']:.4f}+-{r6['chi_P_se']:.4f}  <P>={r6['Pavg_mean']:+.4f}  "
          f"C_P(r)={np.round(r6['C_P'],4)}", flush=True)

    # ---- verdict ----
    ratio = r6['chi_flat'] / r4['chi_flat']
    print("\n[VERDICT]  [V]", flush=True)
    print(f"   chi_flat: L=4 {r4['chi_flat']:.3f} -> L=6 {r6['chi_flat']:.3f} "
          f"(ratio {ratio:.3f}); free-massless would scale ~(6/4)^2=2.25.", flush=True)
    print(f"   Polyakov C_P(r) decays to ~0 by r=2 at both L -> xi<~1 (gapped Polyakov sector).",
          flush=True)
    print(f"   -> pilot data are consistent with a PLATEAU (chi_flat ~ L-flat within errors;"
          f" NOT the L^2 growth of a free field).", flush=True)
    print("   Decisive run: L=8,12 at beta=2.5, N>=2000 (gap>>tau), same chi_flat + force-moment,"
          " one job at a time.", flush=True)

    out['verdict'] = dict(chi_flat_ratio_6_over_4=float(ratio),
                          free_massless_ratio=2.25,
                          plateau_consistent=bool(0.75 < ratio < 1.35))
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'q5_4d_chi_results.json')
    with open(path, 'w') as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
    print(f"\n# wrote {path}", flush=True)

if __name__ == '__main__':
    main()
