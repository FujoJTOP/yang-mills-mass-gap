#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gevp_4d.py -- fill the `[O]` data gap of KT-RHO2 (S7-0): the 4D SU(2) 0^{++} correlator
MATRIX  C_ij(t)  for 3 operators (i,j in {o0=plaquette, o2=APE2-plaquette, o4=APE4-plaquette}),
then the 2-state GEVP  C(t) v = lam(t,t0) C(t0) v  with t0 = 1,2,
    m_n(t,t0) = -ln lam_n /(t-t0),   rho^2(t,t0) = (m2/m1)^2 ,  Q = (rho^2-1)/3 .
Plus the ACID CONTROL of the KT criterion: the SAME code on
    (i) a free massless Gaussian gauge field (g=0 limit, exact Gaussian ensemble),
    (ii) a beta=0 (Haar) field.

WHAT IS ADDED, WHAT IS NOT.
  * The sampler is NOT touched: make_stream / plaq_operator_t / plaq are imported verbatim
    from m1_4d_run.py (checkerboard heat bath, identical protocol therm/gap, cold start).
  * The only addition is the accumulation and the on-disk persistence of the CROSS
    correlators C_ij(t)  (m1_4d_run computed C0/C2 diagonals and discarded them, and stored
    only C_P).
  * No change of physics, of the ensemble, or of the operator definitions.

Crash safety.  Every `ckpt_every` stored configurations the operator array (N,nops,L), the
gauge field U and a state JSON are written; a re-launch with the same (L,N,protocol) resumes
at n_done with the saved field (in-place U[:] = saved).

Usage (py -3.12 -I, numpy/scipy only):
    py -3.12 gevp_4d.py measure 4 400          # smoke
    py -3.12 gevp_4d.py measure 8 4000         # production
    py -3.12 gevp_4d.py analyze 4              # -> gevp_4d_L4.json
    py -3.12 gevp_4d.py free 8 2000            # free-field control -> gevp_4d_free_L8.json
    py -3.12 gevp_4d.py haar 8 400             # beta=0 control
Labels [V]; no QED.  Clay 0%.
"""
import sys, os, json, time
import numpy as np

sys.path.insert(0, r'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, r'D:/Dev/math/research/unconditional/w7')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import w4_blocked as W          # noqa: E402
import m1_4d_run as M           # noqa: E402  (sampler + operators, UNCHANGED)

BETA = 2.5
ALPHA = 0.5
STEPS = (2, 4)                  # APE steps for o2, o4  (KT-RHO2 S7-0 item 1)
TAGS = ['o0_plaq', 'o2_ape2', 'o4_ape4']
PROTO = 'therm=800,gap=10,cold,beta=2.5,alpha=0.5,APE2+APE4,SU(2)periodic'


def paths(L):
    return dict(ops=os.path.join(HERE, f'gevp_4d_L{L}_ops.npy'),
                field=os.path.join(HERE, f'gevp_4d_L{L}_field.npy'),
                plq=os.path.join(HERE, f'gevp_4d_L{L}_plq.npy'),
                state=os.path.join(HERE, f'gevp_4d_L{L}_state.json'))


# ==========================================================================
# 1. measurement   (sampler untouched -- imported)
# ==========================================================================
def observe(U):
    """The 3 A1++ zero-momentum operators of this configuration, shape (nops, L)."""
    out = [M.plaq_operator_t(U)]
    for s in STEPS:
        out.append(M.plaq_operator_t(W.ape_all2(U, ALPHA, s)))
    return np.array(out)


def measure(L, N, therm=800, gap=10, seed=101, ckpt_every=250, resume=True, quiet=False):
    P = paths(L)
    proto = f'{PROTO},therm={therm},gap={gap}'
    nops = len(TAGS)
    U, sweep = M.make_stream(L, BETA, seed, cold=True)
    A = np.zeros((N, nops, L))
    plq = np.zeros(N)
    n_done = 0
    wall0 = time.time()
    if resume and os.path.exists(P['state']) and os.path.exists(P['ops']) \
            and os.path.exists(P['field']):
        try:
            st = json.load(open(P['state']))
            if st.get('L') == L and st.get('N') == N and st.get('protocol') == proto:
                Ud = np.load(P['field'])
                Ad = np.load(P['ops'])
                if Ud.shape == U.shape and Ad.shape[1:] == (nops, L):
                    U[:] = Ud
                    n_done = min(Ad.shape[0], N)
                    A[:n_done] = Ad[:n_done]
                    if os.path.exists(P['plq']):
                        qd = np.load(P['plq'])
                        plq[:min(len(qd), N)] = qd[:min(len(qd), N)]
                    print(f'# resume L={L} at n_done={n_done}', flush=True)
        except Exception as e:
            print(f'# resume failed ({e!r}); fresh chain', flush=True)
            n_done = 0
    if n_done == 0:
        t0 = time.time()
        for _ in range(therm):
            sweep()
        if not quiet:
            print(f'# L={L} therm {therm} sweeps in {time.time()-t0:.1f}s', flush=True)
    tstart = time.time()
    for i in range(n_done, N):
        for _ in range(gap):
            sweep()
        A[i] = observe(U)
        plq[i] = M.plaq(U, BETA)
        if ((i + 1) % ckpt_every == 0) or (i + 1 == N):
            np.save(P['ops'], A[:i + 1])
            np.save(P['field'], U)
            np.save(P['plq'], plq[:i + 1])
            nd = i + 1 - n_done
            json.dump(dict(L=L, N=N, n_done=int(i + 1), protocol=proto, beta=BETA,
                           alpha=ALPHA, steps=list(STEPS), tags=TAGS, seed=seed,
                           therm=therm, gap=gap, wall_s=time.time() - wall0,
                           s_per_conf=(time.time() - tstart) / max(1, nd)),
                      open(P['state'], 'w'), indent=1)
            if not quiet:
                print(f'  L={L} n={i+1}/{N}  plaq={plq[i]:.5f}  '
                      f'({(time.time()-tstart)/max(1,nd):.3f} s/conf)', flush=True)
    return A, plq, dict(L=L, N=int(N), n_done=int(N), therm=therm, gap=gap, seed=seed,
                        wall_s=time.time() - wall0)


# ==========================================================================
# 2. correlator matrix + GEVP + bootstrap
# ==========================================================================
def per_config_outer(ops):
    """S[n,i,j,t] = sum_{t0} o_i(t0+t) o_j(t0);  abar[n,i] = mean_t o_i."""
    N, no, L = ops.shape
    S = np.empty((N, no, no, L))
    for t in range(L):
        S[:, :, :, t] = np.einsum('nil,njl->nij', ops, np.roll(ops, -t, axis=2))
    return S, ops.mean(axis=2)


def C_from_w(c, S, abar, L):
    """Connected C_ij(t) from config-multiplicity weights c (sum c = N) -- exact.

    Uses the reflection identity C_ij(t) = C_ji(L-t) and returns the symmetrised
    C_ij(t) <- (C_ij(t) + C_ji(L-t))/2 (same expectation, halves the noise coming from the
    time-translation average, which at L=4 rests on 4 shifts only).
    """
    N, no = abar.shape
    Nt = c.sum()
    num = np.tensordot(c, S.reshape(N, -1), axes=(0, 0)).reshape(no, no, L)
    om = (c @ abar) / Nt
    C = num / (Nt * L) - om[:, None, None] * om[None, :, None]
    # reflection: C_ij(t) = C_ji((-t) mod L) = C_ji(L-t)
    D = C.transpose(1, 0, 2)
    refl = np.roll(D[:, :, ::-1], 1, axis=2)      # refl[i,j,t] = C[j,i,(-t)%L]
    return 0.5 * (C + refl)


def gevp(Amat, Bmat, keep_min=2):
    """Eigenvalues lam of C(t) v = lam C(t0) v (symmetric whitening of C(t0))."""
    A = 0.5 * (Amat + Amat.T)
    B = 0.5 * (Bmat + Bmat.T)
    w, V = np.linalg.eigh(B)
    keep = w > 1e-10 * max(abs(w).max(), 1e-300)
    if keep.sum() < keep_min:
        return None, int(keep.sum()), float('nan')
    Wv = V[:, keep] / np.sqrt(w[keep])
    Mx = Wv.T @ A @ Wv
    Mx = 0.5 * (Mx + Mx.T)
    lam = np.linalg.eigh(Mx)[0][::-1]
    cond = float(abs(w[keep]).max() / abs(w[keep]).min())
    return lam, int(keep.sum()), cond


SUBSETS = {'2op_o0o2': (0, 1), '2op_o2o4': (1, 2), '3op_o0o2o4': (0, 1, 2)}


def table_for(Cm, L, idx, label):
    """GEVP table for the operator subset `idx`.

    lam_n(t,t0) solves C(t) v = lam C(t0) v; m_log,n = -ln lam_n /(t-t0)  (this is the
    quantity the KT criterion uses through rho^2 = (m2/m1)^2).  The periodic time extent
    contaminates the plain log at t ~ L/2, so the cosh-corrected variant
        m_cosh,n(t) = arccosh( (lam_n(t+1)+lam_n(t-1)) / (2 lam_n(t)) )
    is reported alongside (needs t+1 < L).
    """
    sub = Cm[np.ix_(list(idx), list(idx), list(range(Cm.shape[2])))]
    rows = []
    for t0 in (1, 2):
        lam_t = {}
        for t in range(t0, L):
            lam_t[t] = gevp(sub[:, :, t], sub[:, :, t0])
        for t in range(t0 + 1, L // 2 + 1):
            lam, rank, cond = lam_t[t]
            row = dict(subset=label, ops=list(idx), t0=t0, t=t, dt=t - t0,
                       rank=rank, cond=cond, ok=False)
            if lam is None or len(lam) < 2 or lam[0] <= 0 or lam[1] <= 0:
                rows.append(row)
                continue
            m1 = float(-np.log(lam[0]) / (t - t0))
            m2 = float(-np.log(lam[1]) / (t - t0))
            rho2 = float((m2 / m1) ** 2) if m1 != 0 else float('nan')
            row.update(ok=True, lam=[float(lam[0]), float(lam[1])], m1=m1, m2=m2,
                       rho2=rho2, Q=float((rho2 - 1.0) / 3.0))
            # cosh-corrected masses
            lo = lam_t.get(t - 1)
            hi = lam_t.get(t + 1)
            if lo is not None and hi is not None and lo[0] is not None and hi[0] is not None:
                mc = []
                for n in range(2):
                    arg = (lo[0][n] + hi[0][n]) / (2.0 * lam[n])
                    mc.append(float(np.arccosh(arg)) if arg >= 1.0 else float('nan'))
                if mc[0] == mc[0] and mc[0] != 0 and mc[1] == mc[1]:
                    row['m1_cosh'] = mc[0]
                    row['m2_cosh'] = mc[1]
                    row['rho2_cosh'] = float((mc[1] / mc[0]) ** 2)
                    row['Q_cosh'] = float((row['rho2_cosh'] - 1.0) / 3.0)
            rows.append(row)
    return rows


def _boot_stats(row, arr):
    """attach bootstrap mean/error/percentiles to a table row."""
    if row is None or not row.get('ok'):
        return
    arr = arr[np.isfinite(arr)]
    if arr.size < 3:
        return
    row['rho2_se'] = float(arr.std(ddof=1))
    row['rho2_lo'] = float(np.percentile(arr, 16))
    row['rho2_hi'] = float(np.percentile(arr, 84))
    row['Q_se'] = float(arr.std(ddof=1) / 3.0)
    row['frac_gt4'] = float(np.mean(arr > 4.0))
    row['frac_gt225'] = float(np.mean(arr > 2.25))
    row['frac_positive'] = float(np.mean(arr > 0.0))
    row['n_boot_ok'] = int(arr.size)


def analyze(L, ops=None, plq=None, meta=None, boots=400, tag=''):
    if ops is None:
        ops = np.load(paths(L)['ops'])
    N, no, Lz = ops.shape
    assert Lz == L, (Lz, L)
    S, abar = per_config_outer(ops)
    C0 = C_from_w(np.ones(N), S, abar, L)

    d = np.sqrt(np.maximum(np.diag(C0[:, :, 0]), 1e-300))
    G = C0[:, :, 0] / np.outer(d, d)
    evG = np.sort(np.linalg.eigvalsh(0.5 * (G + G.T)))[::-1]

    # ---- bootstrap (nonparametric over configurations) ----
    rng = np.random.default_rng(12345)
    acc = {}
    for _ in range(boots):
        c = np.bincount(rng.integers(0, N, N), minlength=N).astype(float)
        Cb = C_from_w(c, S, abar, L)
        for lab, idx in SUBSETS.items():
            for row in table_for(Cb, L, idx, lab):
                if row['ok']:
                    acc.setdefault((lab, row['t0'], row['t']), []).append(row['rho2'])
                    if 'rho2_cosh' in row:
                        acc.setdefault((lab, row['t0'], row['t'], 'cosh'), []).append(
                            row['rho2_cosh'])

    tables = {}
    for lab, idx in SUBSETS.items():
        rows = table_for(C0, L, idx, lab)
        for row in rows:
            _boot_stats(row, np.array(acc.get((lab, row['t0'], row['t']), []), float))
            ch = acc.get((lab, row['t0'], row['t'], 'cosh'))
            if ch is not None and 'rho2_cosh' in row:
                a = np.array(ch, float)
                a = a[np.isfinite(a)]
                if a.size > 2:
                    row['rho2_cosh_se'] = float(a.std(ddof=1))
                    row['rho2_cosh_lo'] = float(np.percentile(a, 16))
                    row['rho2_cosh_hi'] = float(np.percentile(a, 84))
        tables[lab] = rows

    diag_meff = [[(float(np.log(C0[i, i, t] / C0[i, i, t + 1]))
                   if C0[i, i, t] > 0 and C0[i, i, t + 1] > 0 else float('nan'))
                  for t in range(1, max(2, L // 2))] for i in range(no)]
    diag_meff_cosh = []
    for i in range(no):
        row = []
        for t in range(1, max(2, L // 2)):
            a_, b_, c_ = C0[i, i, t - 1], C0[i, i, t], C0[i, i, t + 1]
            arg = (a_ + c_) / (2.0 * b_) if b_ > 0 else -1.0
            row.append(float(np.arccosh(arg)) if arg >= 1.0 else float('nan'))
        diag_meff_cosh.append(row)

    md = dict(L=L, N=int(N), nops=int(no), tags=TAGS, beta=BETA, alpha=ALPHA,
              steps=list(STEPS), protocol=PROTO, script='gevp_4d.py',
              bootstrap=dict(kind='nonparametric over configs', B=boots, seed=12345),
              readout='rho^2=(m2/m1)^2, m_n=-ln lam_n/(t-t0), Q=(rho^2-1)/3; '
                      'cosh-corrected variants also reported')
    md.update(meta or {})
    md.update(L=L, N=int(N), nops=int(no))
    out = dict(meta=md,
               C_t=[[[float(C0[i, j, t]) for t in range(L)] for j in range(no)]
                    for i in range(no)],
               C_t_scaled=[[[float(C0[i, j, t] /
                                   np.sqrt(max(C0[i, i, 0] * C0[j, j, 0], 1e-300)))
                             for t in range(L)] for j in range(no)] for i in range(no)],
               diag_meff=diag_meff, diag_meff_cosh=diag_meff_cosh,
               Gram_eigs=[float(x) for x in evG], gram_condition=float(evG[0] / evG[-1]),
               tables=tables)
    if plq is not None and len(plq) > 1:
        out['plq'] = dict(mean=float(np.mean(plq)), sd=float(np.std(plq, ddof=1)),
                          se=float(np.std(plq, ddof=1) / np.sqrt(len(plq))), n=int(len(plq)))
    if tag:
        out['meta']['tag'] = tag
    return out


# ==========================================================================
# 3. controls: free massless Gaussian gauge field, beta=0 Haar
# ==========================================================================
def free_su2_config(L, rng, g=0.30, mass=0.0):
    """One Gaussian free gauge-field configuration (Landau-transverse, k=0 mode dropped).

    Gaussian measure with (1/(4g^2)) sum tr F^2  =>  propagator
        g^2 (delta_munu - khat_mu khat_nu / khat^2) / (khat^2 + m^2) ,
    khat_mu = 2 sin(p_mu/2), p_mu = 2 pi n_mu / L.  Links  U_mu = exp(i A^a_mu sigma^a/2)
    are exactly SU(2) since A^a_mu is real (and real by FFT Hermiticity).
    """
    ks = np.stack(np.meshgrid(*[np.arange(L)] * 4, indexing='ij'), axis=-1).astype(float)
    kh = 2.0 * np.sin(np.pi * ks / L)
    k2 = (kh ** 2).sum(-1)
    k2s = np.where(k2 > 0, k2, 1.0)
    w = rng.standard_normal((3,) + (L,) * 4 + (4,))
    F = np.fft.fftn(w, axes=(1, 2, 3, 4)) / L ** 2          # E|F|^2 = N_modes normalised
    dot = np.einsum('a...m,...m->a...', F, kh)
    At = F - kh[None] * (dot / k2s)[..., None]              # transverse projector
    At = At * (g / np.sqrt(k2s + mass ** 2))[None][..., None]
    At[:, 0, 0, 0, 0, :] = 0.0                              # remove the constant mode
    # A[a][x0,x1,x2,x3,mu] , shape (3, L,L,L,L, 4)
    A = np.real(np.fft.ifftn(At, axes=(1, 2, 3, 4))) * L ** 4
    phi = np.sqrt((A ** 2).sum(0))                          # (L,L,L,L,4)
    th = 0.5 * phi
    inv = np.where(phi > 1e-300, 1.0 / np.maximum(phi, 1e-300), 0.0)
    nx, ny, nz = A[0] * inv, A[1] * inv, A[2] * inv
    c, s = np.cos(th), np.sin(th)
    U = np.empty((4, L, L, L, L, 2, 2), dtype=np.complex128)

    def to_u(x):            # components are (L,L,L,L,4) ; U wants (4,L,L,L,L)
        return np.moveaxis(x, -1, 0)
    U[..., 0, 0] = to_u(c + 1j * nz * s)
    U[..., 0, 1] = to_u((1j * nx + ny) * s)
    U[..., 1, 0] = to_u((1j * nx - ny) * s)
    U[..., 1, 1] = to_u(c - 1j * nz * s)
    return U


def measure_free(L, N, g=0.30, seed=7, ckpt_every=500):
    rng = np.random.default_rng(seed)
    A = np.zeros((N, len(TAGS), L))
    plq = np.zeros(N)
    t0 = time.time()
    for i in range(N):
        U = free_su2_config(L, rng, g=g)
        A[i] = observe(U)
        plq[i] = M.plaq(U, BETA)
        if (i + 1) % ckpt_every == 0:
            print(f'  free L={L} n={i+1}/{N} ({(time.time()-t0)/(i+1):.4f} s/conf)', flush=True)
            np.save(os.path.join(HERE, f'gevp_4d_free_L{L}_ops.npy'), A[:i + 1])
    np.save(os.path.join(HERE, f'gevp_4d_free_L{L}_ops.npy'), A)
    return A, plq, dict(L=L, N=int(N), g=g, wall_s=time.time() - t0)


def measure_haar(L, N, seed=11):
    """beta=0: links i.i.d. Haar; plaquettes at different t share no link => C(t>0)=0."""
    U, sweep = M.make_stream(L, 0.0, seed, cold=True)
    A = np.zeros((N, len(TAGS), L))
    for i in range(N):
        for _ in range(2):
            sweep()
        A[i] = observe(U)
    np.save(os.path.join(HERE, f'gevp_4d_haar_L{L}_ops.npy'), A)
    return A, np.zeros(N), dict(L=L, N=int(N), beta=0.0)


# ==========================================================================
def save_json(obj, name):
    def d(o):
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)
    p = os.path.join(HERE, name)
    json.dump(obj, open(p, 'w'), indent=1, default=d)
    return p


def show(res):
    m = res['meta']
    pl = res.get('plq', {})
    print(f"\n# L={m['L']} N={m['N']} beta={m['beta']}  "
          f"plaq={pl.get('mean', float('nan')):.5f}  tag={m.get('tag','YM')}")
    print(f"  Gram eigs of normalised C(0): {np.round(res['Gram_eigs'], 5)}"
          f"   cond={res['gram_condition']:.1f}")
    for i, tg in enumerate(TAGS):
        print(f'  m_eff diag {tg:9s} log={np.round(res["diag_meff"][i], 4)} '
              f'cosh={np.round(res["diag_meff_cosh"][i], 4)}')
    for lab, rows in res['tables'].items():
        print(f"  -- {lab} --")
        print(f"  {'t0':>2} {'t':>2} {'rk':>3} {'lam1':>9} {'lam2':>9} {'m1':>7} {'m2':>7} "
              f"{'rho2':>9} {'+-':>8} {'Q':>8} {'rho2_cosh':>10}")
        for r in rows:
            if not r['ok']:
                print(f"  {r['t0']:2d} {r['t']:2d} {r['rank']:3d}   -- rank deficient --")
                continue
            rc = r.get('rho2_cosh', float('nan'))
            print(f"  {r['t0']:2d} {r['t']:2d} {r['rank']:3d} {r['lam'][0]:9.5f} "
                  f"{r['lam'][1]:9.5f} {r['m1']:7.4f} {r['m2']:7.4f} {r['rho2']:9.4f} "
                  f"{r.get('rho2_se', float('nan')):8.4f} {r['Q']:8.4f} {rc:10.4f}")


if __name__ == '__main__':
    a = sys.argv[1:]
    cmd = a[0] if a else 'help'
    if cmd == 'measure':
        L = int(a[1]); N = int(a[2])
        therm = int(a[3]) if len(a) > 3 else 800
        gap = int(a[4]) if len(a) > 4 else 10
        A, plq, meta = measure(L, N, therm=therm, gap=gap,
                               seed={4: 101, 6: 202, 8: 303, 16: 505}.get(L, 100 + L))
        res = analyze(L, ops=A, plq=plq, meta=meta)
        show(res)
        print('# ->', save_json(res, f'gevp_4d_L{L}.json'))
    elif cmd == 'analyze':
        L = int(a[1])
        pl = paths(L)['plq']
        plq = np.load(pl) if os.path.exists(pl) else None
        st = paths(L)['state']
        meta = json.load(open(st)) if os.path.exists(st) else {}
        res = analyze(L, plq=plq, meta=meta)
        show(res)
        print('# ->', save_json(res, f'gevp_4d_L{L}.json'))
    elif cmd == 'free':
        L = int(a[1]); N = int(a[2])
        g = float(a[3]) if len(a) > 3 else 0.30
        A, plq, meta = measure_free(L, N, g=g)
        res = analyze(L, ops=A, plq=plq, meta=meta, tag=f'free_massless_g={g}')
        show(res)
        print('# ->', save_json(res, f'gevp_4d_free_L{L}.json'))
    elif cmd == 'haar':
        L = int(a[1]); N = int(a[2])
        A, plq, meta = measure_haar(L, N)
        res = analyze(L, ops=A, plq=None, meta=meta, tag='beta=0_Haar')
        show(res)
        print('# ->', save_json(res, f'gevp_4d_haar_L{L}.json'))
    else:
        print(__doc__)
