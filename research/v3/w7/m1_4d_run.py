#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
m1_4d_run.py -- THE decisive 4D SU(2) numerical run for the Yang-Mills campaign.

Target (H1 / carrier, Q5-4D-CHI + SM6 + ch12-numerics): at beta=2.5, periodic 4D SU(2):
  (A) chi = Var(a)  (k=0 susceptibility, a = N^{-1/2} sum_l log(U_l), colour-averaged),
      measured in the FLAT chart (chi_flat, continues the L=4,6 pilot) AND in LANDAU gauge
      (chi_landau, the standard zero-momentum gluon susceptibility) with Gribov-copy control;
  (B) the gauge-invariant 0^{++} correlator  C(t)=<O(t)O(0)>_c,  O=space-summed plaquette
      (unsmeared + APE-smeared),  fit m1 and the residue/f1^2 proxy;
  (C) the gauge-invariant Polyakov spatial correlator C_P(r) -> xi (N-ality-1 cross-check).

Volumes L=4,6,8 (L=12 optional), N>=2000 configs, gap>>tau, ONE sequential process
(one heavy job at a time; charter s4 -- never parallelise SU(N) L>=6).  Results are
written INCREMENTALLY to m1_4d_run_results.json after every completed stage, so partial
work is never lost.  py -3.12 -I  (numpy/scipy).  NO third-party sources.

Labels [T]/[C]/[NG]/[O]/[V]; no QED.
"""
import sys, os, json, time
import numpy as np
from scipy.optimize import least_squares

CPB = 'D:/Dev/math/research/lattice-cpb'
W7U = 'D:/Dev/math/research/unconditional/w7'
HERE = os.path.dirname(os.path.abspath(__file__))
for p in (CPB, W7U):
    sys.path.insert(0, p)
import w4_blocked as W
import gpp_core as C

RESULT_PATH = os.path.join(HERE, 'm1_4d_run_results.json')
SIG = W.SIG
ID2 = W.ID2


# --------------------------------------------------------------------------
# streaming SU(2) heat bath (identical protocol to the Q5 pilot)
# --------------------------------------------------------------------------
def make_stream(L, beta, seed, cold=False):
    rng = np.random.default_rng(seed)
    if cold:
        U = np.zeros((4, L, L, L, L, 2, 2), dtype=np.complex128)
        U[..., 0, 0] = 1.0
        U[..., 1, 1] = 1.0
    else:
        U = W.random_su2(rng, (4, L, L, L, L))
    ns = L ** 4
    par = (np.indices((L, L, L, L)).sum(0)).reshape(-1) % 2
    idxp = [np.nonzero(par == c)[0] for c in (0, 1)]
    kappa = beta / 2.0

    def sweep():
        Uf = U.reshape(4, ns, 2, 2)
        for c in (0, 1):
            ids = idxp[c]
            for mu in range(4):
                S = W.staple_sum2(U, mu).reshape(ns, 2, 2)[ids]
                A = Uf[mu, ids]
                Uf[mu, ids] = W.su2_heatbath(kappa, A @ S, 0, 1, rng) @ A
    return U, sweep


# --------------------------------------------------------------------------
# estimators
# --------------------------------------------------------------------------
def k0_vector(U):
    """A = N^{-1/2} sum_l xi_l in R^3, xi=log U (flat-chart k=0 mode), N=4 L^4."""
    x = C.su2_log_x(U)
    Nl = 4 * U.shape[1] ** 4
    return x.reshape(Nl, 3).sum(0) / np.sqrt(Nl)


def chi_flat(U):
    """chi_flat = (1/3)<|A|^2>, A=N^{-1/2} sum_l xi_l, xi=log U (flat-chart k=0 mode)."""
    A = k0_vector(U)
    return float(A @ A / 3.0)


def polyakov_field(U):
    L = U.shape[1]
    A = np.broadcast_to(np.eye(2, dtype=complex), (L, L, L, L, 2, 2)).copy()
    for t in range(L):
        A = A @ np.roll(U[0], -t, axis=0)
    return np.real(np.trace(A[0], axis1=-2, axis2=-1)) / 2.0


def plaq_operator_t(U):
    """O(t) = sum_{spatial i<j} sum_space (1 - 1/2 Re tr U_ij), per timeslice t."""
    L = U.shape[1]
    op = np.zeros(L)
    for i in range(1, 4):
        for j in range(i + 1, 4):
            P = U[i] @ np.roll(U[j], -1, axis=i) @ W.adj(np.roll(U[i], -1, axis=j)) @ W.adj(U[j])
            pl = 1.0 - np.real(np.trace(P, axis1=-2, axis2=-1)) / 2.0
            op += pl.sum(axis=(1, 2, 3))
    return op


def plaq(U, beta):
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


def jk_ms_var(A, nblk=20):
    """Delete-block jackknife of V = (1/3)(<|A|^2> - |<A>|^2), A (N,3) real.
    Returns (V, SE, mean_contrib=|mean A|^2/3)."""
    A = np.asarray(A, float); n = A.shape[0]
    nblk = int(min(nblk, max(4, n // 8)))
    ed = np.linspace(0, n, nblk + 1).astype(int)

    def Vof(Ab):
        return float((Ab ** 2).sum(1).mean() / 3.0 - (Ab.mean(0) @ Ab.mean(0)) / 3.0)

    base = Vof(A)
    lv = []
    for b in range(nblk):
        idx = np.ones(n, bool); idx[ed[b]:ed[b + 1]] = False
        lv.append(Vof(A[idx]))
    lv = np.array(lv)
    se = float(np.sqrt((nblk - 1) / nblk * np.sum((lv - lv.mean()) ** 2)))
    return base, se, float(np.mean(A, 0) @ np.mean(A, 0) / 3.0)


def jk_mean_se(x, nblk=20):
    x = np.asarray(x, float); n = x.size
    nblk = int(min(nblk, max(4, n // 8)))
    ed = np.linspace(0, n, nblk + 1).astype(int)
    m = x.mean()
    lv = np.array([np.delete(x, np.arange(ed[b], ed[b + 1])).mean() for b in range(nblk)])
    se = np.sqrt((nblk - 1) / nblk * np.sum((lv - lv.mean()) ** 2))
    return float(m), float(se)


# --------------------------------------------------------------------------
# Landau gauge fixing -- red-black, exact SU(2) local optimizer + Gribov control
# --------------------------------------------------------------------------
def project_su2(D):
    """SU(2) g maximizing Re tr(g D): D (...,2,2) -> g.  Exact for SU(2)."""
    d0 = 0.5 * np.trace(D, axis1=-2, axis2=-1).real
    d = 0.5 * np.einsum('aij,...ji->...a', SIG, D)
    v = np.concatenate([d0[..., None], -d.imag], axis=-1)
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    n = np.where(n < 1e-14, 1.0, n)
    q = v / n
    return q[..., 0, None, None] * ID2 + 1j * np.einsum('...a,aij->...ij', q[..., 1:], SIG)


def gauge_transform(U, G):
    out = np.empty_like(U)
    for mu in range(4):
        out[mu] = G @ U[mu] @ W.adj(np.roll(G, -1, axis=mu))
    return out


def landau_functional(V):
    """F = sum_{x,mu} Re tr V_mu(x) (maximised by Landau gauge)."""
    return float(sum(np.trace(V[mu], axis1=-2, axis2=-1).real.sum() for mu in range(4)))


def landau_fix(U, tol=1e-7, maxit=400):
    """Red-black over-relaxation-free Landau gauge fixing. Returns (V, F, iters)."""
    L = U.shape[1]
    V = [U[mu].copy() for mu in range(4)]
    shape = (L, L, L, L)
    par = (np.indices(shape).sum(0)) % 2
    eye = np.broadcast_to(ID2, shape + (2, 2))
    its = maxit
    for it in range(maxit):
        for c in (0, 1):
            D = np.zeros(shape + (2, 2), complex)
            for mu in range(4):
                D += V[mu] + W.adj(np.roll(V[mu], 1, axis=mu))
            Gopt = project_su2(D)
            mask = (par == c)
            G = np.where(mask[..., None, None], Gopt, eye)
            for mu in range(4):
                V[mu] = G @ V[mu] @ W.adj(np.roll(G, -1, axis=mu))
        D = np.zeros(shape + (2, 2), complex)
        for mu in range(4):
            D += V[mu] + W.adj(np.roll(V[mu], 1, axis=mu))
        g0 = np.trace(project_su2(D), axis1=-2, axis2=-1).real / 2
        dev = float(np.mean(np.arccos(np.clip(g0, -1, 1))))
        if it > 2 and dev < tol:
            its = it + 1
            break
    Vout = np.stack(V)
    return Vout, landau_functional(Vout), its


def random_su2_field(rng, shape):
    return W.random_su2(rng, shape)


# --------------------------------------------------------------------------
# 0^{++} correlator fit: C(t)=A [ exp(-m t) + exp(-m (L-t)) ]
# --------------------------------------------------------------------------
def eff_mass(Cc, L):
    """m_eff(t) = log(C(t)/C(t+1)) for t=1..L//2-1 (model-independent sanity)."""
    out = []
    for t in range(1, L // 2):
        if Cc[t] > 0 and Cc[t + 1] > 0:
            out.append(float(np.log(Cc[t] / Cc[t + 1])))
        else:
            out.append(float('nan'))
    return out


def fit_corr_length(C_P, L):
    """Fit C_P(r) = C_P(0) exp(-r/xi) over r=1..L//2 (positive). Returns xi."""
    r = np.arange(1, L // 2 + 1)
    y = np.asarray(C_P, float)[r]
    m = y > 0
    if m.sum() < 2:
        return float('nan')
    rr = r[m]; yy = np.log(y[m])
    A = np.vstack([np.ones_like(rr), -rr]).T      # log C = a - r/xi
    coef, *_ = np.linalg.lstsq(A, yy, rcond=None)
    xi = 1.0 / coef[1] if coef[1] > 1e-9 else float('inf')
    return float(xi)


def fit_glueball(Cc, L, tmax=None):
    """Fit connected C(t), t=1..L//2 (drop contact t=0), to A(exp(-m t)+exp(-m(L-t)))."""
    tt = np.arange(1, L // 2 + 1)
    if tmax is not None:
        tt = tt[tt <= tmax]
    y = np.asarray(Cc)[tt]
    if np.any(y <= 0):
        y = np.maximum(y, 1e-30)
    def resid(p):
        A, m = np.exp(p[0]), np.exp(p[1])
        model = A * (np.exp(-m * tt) + np.exp(-m * (L - tt)))
        return (np.log(model) - np.log(y))
    try:
        sol = least_squares(resid, x0=[np.log(max(y[0], 1e-6)), np.log(0.6)],
                            bounds=([-60, np.log(0.03)], [60, np.log(3.0)]), max_nfev=4000)
        A, m = float(np.exp(sol.x[0])), float(np.exp(sol.x[1]))
        rms = float(np.sqrt(np.mean(sol.fun ** 2)))
        return dict(m1=m, A=A, f1sq_proxy=A, logfit_rms=rms, ok=True)
    except Exception as e:
        return dict(m1=float('nan'), A=float('nan'), f1sq_proxy=float('nan'),
                    error=repr(e), ok=False)


# --------------------------------------------------------------------------
# one volume
# --------------------------------------------------------------------------
def run_volume(L, beta, N, therm, gap, seed, n_gf, n_gribov, n_copy, rng_seed):
    U, sweep = make_stream(L, beta, seed, cold=True)
    t0 = time.time()
    for _ in range(therm):
        sweep()
    therm_s = time.time() - t0

    chF = np.empty(N); plq = np.empty(N); Pavg = np.empty(N)
    Cp = np.zeros(L)
    sum_o0 = np.zeros(L); sum_o2 = np.zeros(L)
    C_raw0 = np.zeros(L); C_raw2 = np.zeros(L)
    n_corr = 0
    chi_landau = []
    landau_avec = []
    Aflat = np.empty((N, 3))
    gribov = []
    P0 = None

    gf_stride = max(1, N // max(1, n_gf))
    gribov_ids = set(int(round(k * (N - 1) / max(1, n_gribov - 1))) for k in range(n_gribov)) \
        if n_gribov > 0 else set()
    rng_g = np.random.default_rng(rng_seed)

    def accumulate_corr(o, Craw, sum_o):
        Lc = len(o)
        for t0i in range(Lc):
            oc = np.roll(o, -t0i)
            Craw += o * oc            # sum_t0 o[(t0+dt)%L] o[t0]
        sum_o += o
        return Craw, sum_o

    for i in range(N):
        for _ in range(gap):
            sweep()
        Av = k0_vector(U)
        Aflat[i] = Av
        chF[i] = float(Av @ Av / 3.0)
        plq[i] = plaq(U, beta)
        P = polyakov_field(U)
        Pavg[i] = P.mean()
        Cp += np.array([float(np.mean(P * np.roll(P, -r, axis=0))) for r in range(L)])
        # glueball operators
        o0 = plaq_operator_t(U)
        U2 = W.ape_all2(U, alpha=0.5, nsteps=2)
        o2 = plaq_operator_t(U2)
        C_raw0, sum_o0 = accumulate_corr(o0, C_raw0, sum_o0)
        C_raw2, sum_o2 = accumulate_corr(o2, C_raw2, sum_o2)
        n_corr += 1
        if P0 is None:
            P0 = P.copy()
        # Landau gauge fixing (subsample)
        if (i % gf_stride == 0) and (len(chi_landau) < n_gf):
            G = random_su2_field(rng_g, (L, L, L, L))
            Ufix, F, its = landau_fix(gauge_transform(U, G))
            av_fix = k0_vector(Ufix)
            landau_avec.append(av_fix)
            chi_landau.append(float(av_fix @ av_fix / 3.0))
            if i in gribov_ids and len(gribov) < n_gribov:
                vals = []
                bestF = F; bestchi = chi_landau[-1]
                for _c in range(n_copy - 1):
                    G2 = random_su2_field(rng_g, (L, L, L, L))
                    Uf2, F2, _ = landau_fix(gauge_transform(U, G2))
                    ch2 = chi_flat(Uf2)
                    vals.append(ch2)
                    if F2 > bestF:
                        bestF, bestchi = F2, ch2
                allv = [chi_landau[-1]] + vals
                gribov.append(dict(chi_copies=allv, spread=float(np.std(allv)),
                                   spread_rel=float(np.std(allv) / max(1e-9, np.mean(allv))),
                                   best=bestchi))

    Cp /= N
    Pavg_mean = float(Pavg.mean())
    C_P = Cp - Pavg_mean ** 2
    chi_P = float(L ** 3 * Pavg.var())

    # connected glueball correlators
    Om0 = float(sum_o0.mean()) / n_corr if n_corr else 0.0
    Om2 = float(sum_o2.mean()) / n_corr if n_corr else 0.0
    C0 = C_raw0 / n_corr - Om0 ** 2
    C2 = C_raw2 / n_corr - Om2 ** 2
    fit0 = fit_glueball(C0, L)
    fit2 = fit_glueball(C2, L)

    mF, sF = jk_mean_se(chF); mm, _ = jk_mean_se(plq)
    tF = tau_int(chF)
    ms_f, ms_f_se, mean_contrib_f = jk_ms_var(Aflat)
    res = dict(L=L, beta=beta, N=N, therm=therm, gap=gap, seed=seed,
               wall_therm_s=therm_s, wall_total_s=time.time() - t0,
               chi_flat=mF, chi_flat_se=sF, chi_flat_tau=tF, chi_flat_neff=N / tF,
               chi_flat_sd=float(np.std(chF)),
               chi_flat_ms=ms_f, chi_flat_ms_se=ms_f_se, chi_flat_mean_contrib=mean_contrib_f,
               plaq=mm, Pavg_mean=Pavg_mean, chi_P=chi_P,
               C_P=[float(c) for c in C_P], xi_P=fit_corr_length(C_P, L),
               glueball=dict(O_unsmeared=C0.tolist(), O_ape2=C2.tolist(),
                             Om_unsmeared=Om0, Om_ape2=Om2,
                             fit_unsmeared=fit0, fit_ape2=fit2,
                             meff_unsmeared=eff_mass(C0, L), meff_ape2=eff_mass(C2, L)))
    if chi_landau:
        cl = np.array(chi_landau)
        clm, clse = jk_mean_se(cl)
        Lv = np.array(landau_avec)
        lms, lms_se, lmean = jk_ms_var(Lv)
        res['landau'] = dict(n=len(cl), chi_mean=clm, chi_se=clse,
                             chi_sd=float(np.std(cl)),
                             chi_ms=lms, chi_ms_se=lms_se, mean_contrib=lmean,
                             chi_tau=tau_int(cl))
    if gribov:
        spreads = np.array([g['spread'] for g in gribov])
        rel = np.array([g['spread_rel'] for g in gribov])
        best = np.array([g['best'] for g in gribov])
        res['gribov'] = dict(n_copies=n_copy, n_configs=len(gribov),
                             mean_spread=float(spreads.mean()),
                             mean_spread_rel=float(rel.mean()),
                             max_spread_rel=float(rel.max()),
                             best_mean=float(best.mean()),
                             copies=gribov[:24])
    return res


# --------------------------------------------------------------------------
# incremental JSON writer
# --------------------------------------------------------------------------
def save(out):
    with open(RESULT_PATH, 'w') as fh:
        json.dump(out, fh, indent=1,
                  default=lambda o: int(o) if isinstance(o, np.integer)
                  else (float(o) if isinstance(o, (np.floating, np.ndarray)) else str(o)))


# --------------------------------------------------------------------------
def main():
    args = sys.argv[1:]
    quick = 'quick' in args
    vols = [int(a) for a in args if a.isdigit()]
    if not vols:
        vols = [4, 6, 8]
    if quick:
        N, therm, gap = 60, 200, 4
        n_gf, n_gribov, n_copy = 12, 6, 3
    else:
        N, therm, gap = 2000, 800, 10
        n_gf, n_gribov, n_copy = 700, 48, 5

    print(f"# M1 4D run  beta=2.5 SU(2) periodic  vols={vols} N={N} therm={therm} gap={gap} "
          f"n_gf={n_gf} n_gribov={n_gribov} x{n_copy}  one job", flush=True)
    if os.path.exists(RESULT_PATH):
        try:
            out = json.load(open(RESULT_PATH))
        except Exception:
            out = {}
    else:
        out = {}
    out.setdefault('meta', {})
    out['meta'].update(dict(date='2026-11-24', beta=2.5, group='SU(2)',
                            protocol=f'N={N} therm={therm} gap={gap} n_gf={n_gf} '
                                     f'gribov={n_gribov}x{n_copy}',
                            script='m1_4d_run.py'))
    seeds = {4: 101, 6: 202, 8: 303, 12: 404, 16: 505}
    for L in vols:
        print(f"\n=== L={L} : N={N} therm={therm} gap={gap} ===", flush=True)
        r = run_volume(L, 2.5, N, therm, gap, seeds.get(L, 100 + L),
                       n_gf, n_gribov, n_copy, rng_seed=9000 + L)
        out[f'L{L}'] = r
        save(out)
        print(f"  chi_flat={r['chi_flat']:.4f}+-{r['chi_flat_se']:.4f} "
              f"(tau={r['chi_flat_tau']:.2f}, Neff={r['chi_flat_neff']:.0f}, sd={r['chi_flat_sd']:.2f})",
              flush=True)
        print(f"    chi_flat_ms(Var)={r['chi_flat_ms']:.4f}+-{r['chi_flat_ms_se']:.4f} "
              f"mean_contrib={r['chi_flat_mean_contrib']:.3f}", flush=True)
        print(f"  plaq={r['plaq']:.4f}  chi_P={r['chi_P']:.4f}  xi_P={r['xi_P']:.3f}  "
              f"C_P={np.round(r['C_P'],4)}", flush=True)
        if 'landau' in r:
            l = r['landau']
            print(f"  chi_landau=<|A|^2>/3={l['chi_mean']:.4f}+-{l['chi_se']:.4f} (n={l['n']}, tau={l['chi_tau']:.2f})"
                  f"  Var(ms)={l['chi_ms']:.4f}+-{l['chi_ms_se']:.4f} mean_contrib={l['mean_contrib']:.3f}",
                  flush=True)
        if 'gribov' in r:
            g = r['gribov']
            print(f"  Gribov: mean spread={g['mean_spread']:.3f} "
                  f"(rel {g['mean_spread_rel']*100:.1f}%, max {g['max_spread_rel']*100:.1f}%) "
                  f"best mean={g['best_mean']:.3f}", flush=True)
        gb = r['glueball']
        for tag in ('fit_unsmeared', 'fit_ape2'):
            f = gb[tag]
            print(f"  0++ {tag}: m1={f['m1']:.4f} A(f1^2)={f['A']:.3g} rms={f.get('logfit_rms',float('nan')):.3f}",
                  flush=True)
        print(f"  m_eff unsmeared={np.round(gb['meff_unsmeared'],3)} ape2={np.round(gb['meff_ape2'],3)}",
              flush=True)
        print(f"  (wall {r['wall_total_s']:.1f}s)  -> saved {RESULT_PATH}", flush=True)

    # cross-volume ratios
    print("\n=== ratios ===")
    have = [L for L in (4, 6, 8, 12) if f'L{L}' in out]
    for a, b in zip(have, have[1:]):
        ra, rb = out[f'L{a}'], out[f'L{b}']
        ratio = rb['chi_flat'] / ra['chi_flat']
        free = (b / a) ** 2
        print(f"  chi_flat L={a}->{b}: {ra['chi_flat']:.3f}->{rb['chi_flat']:.3f} "
              f"ratio={ratio:.3f}  (free-mss L^2 ratio {free:.2f})")
    out['meta']['done_vols'] = have
    save(out)
    print(f"\n# wrote {RESULT_PATH}")


if __name__ == '__main__':
    main()
