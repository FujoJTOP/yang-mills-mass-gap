#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
hzw_l2_za.py -- HZW-L2: single-point reading of
    Z_a := Var_mu( (1/N_B) sum_{x in B} delta trF2 ),  delta trF2 = a^{-4}(1 - (1/N) Re tr U_plaq)
    varpi := a^4 b^4 Z_a          (a = 1 here, so varpi = b^4 Var(block average of q))
at L=2, SU(2), beta = 2.5.  Also L=4 cross-check point (same code, same beta).

Exact SU(2) heat bath (von Mises-Fisher on S^3), checkerboard parity.
Validation gates:
  V1  beta=0  -> Var(q) = 1/4 exactly (one-plaquette Haar)
  V2  L=4, beta=2.5 -> <P> = 0.6526 (repo value, HT3/HV2-PLAQ)
  V3  L=2, beta=2.5 -> Metropolis vs heat bath agree
  V4  staple identity  Re tr(U S) = sum_{6 plaq containing link} Re tr(P)
"""
import numpy as np
import json, time, os, sys

ID2 = np.eye(2, dtype=complex)
S1 = np.array([[0, 1], [1, 0]], dtype=complex)
S2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
S3 = np.array([[1, 0], [0, -1]], dtype=complex)


def dag(A):
    return A.conj().swapaxes(-1, -2)


def quat_to_su2(q):
    u0, u1, u2, u3 = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return (u0[..., None, None] * ID2
            + 1j * (u1[..., None, None] * S1 + u2[..., None, None] * S2 + u3[..., None, None] * S3))


def random_su2(rng, shape):
    q = rng.standard_normal(shape + (4,))
    q /= np.linalg.norm(q, axis=-1, keepdims=True)
    return quat_to_su2(q)


def parity_sets(L):
    x = np.indices((L,) * 4)
    par = (x.sum(axis=0)) % 2
    return [np.nonzero(par == 0), np.nonzero(par == 1)]


def staple_mu(U, mu, idx):
    """S_mu(x) = sum_{nu != mu} [ U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag
                              + U_nu(x-nu+mu)^dag U_mu(x-nu)^dag U_nu(x-nu) ]   (at sites idx)"""
    S = np.zeros((len(idx[0]), 2, 2), dtype=complex)
    for nu in range(4):
        if nu == mu:
            continue
        Un = U[nu]
        # +nu plaquette based at x
        A = np.roll(Un, -1, axis=mu)[idx]
        B = dag(np.roll(U[mu], -1, axis=nu)[idx])
        C = dag(Un[idx])
        S += A @ B @ C
        # -nu plaquette based at x-nu
        Unm = np.roll(Un, 1, axis=nu)
        Um = np.roll(U[mu], 1, axis=nu)
        A2 = dag(np.roll(Unm, -1, axis=mu)[idx])
        B2 = dag(Um[idx])
        C2 = Unm[idx]
        S += A2 @ B2 @ C2
    return S


def _sample_w0(lam, rng):
    """density propto sqrt(1-x^2) exp(lam x) on [-1,1]  (exact, rejection on inversion)."""
    lam = np.asarray(lam, dtype=np.float64)
    x = np.empty(lam.shape[0])
    todo = np.ones(lam.shape[0], dtype=bool)
    l = np.maximum(lam, 1e-9)
    while todo.any():
        m = todo.sum()
        r = rng.random(m)
        v = -1.0 + np.log1p(r * np.expm1(2.0 * l[todo])) / l[todo]
        v = np.clip(v, -1.0, 1.0)
        acc = rng.random(m) < np.sqrt(np.maximum(1.0 - v * v, 0.0))
        ids = np.nonzero(todo)[0]
        x[ids[acc]] = v[acc]
        todo[ids[acc]] = False
    return x


def heatbath_link(U, mu, idx, beta, rng):
    S = staple_mu(U, mu, idx)
    A = U[mu][idx]
    M = A @ S
    m0 = 0.5 * (M[:, 0, 0] + M[:, 1, 1])
    m1 = 0.5 * (M[:, 0, 1] + M[:, 1, 0])
    m2 = 0.5j * (M[:, 0, 1] - M[:, 1, 0])
    m3 = 0.5 * (M[:, 0, 0] - M[:, 1, 1])
    kappa = beta / 2.0
    a = kappa * np.stack([2 * m0.real, -2 * m1.imag, -2 * m2.imag, -2 * m3.imag], axis=1)
    lam = np.linalg.norm(a, axis=1)
    n = len(lam)
    w0 = _sample_w0(lam, rng)
    g = rng.standard_normal((n, 3))
    g /= np.linalg.norm(g, axis=1, keepdims=True)
    s = np.sqrt(np.maximum(1.0 - w0 * w0, 0.0))
    W = np.concatenate([w0[:, None], s[:, None] * g], axis=1)
    small = lam < 1e-12
    ah = np.where(small[:, None], np.array([1.0, 0.0, 0.0, 0.0]), a / np.maximum(lam, 1e-300)[:, None])
    # quaternion product v = ah * W  (rotation of the aligned sample)
    v0 = ah[:, 0] * W[:, 0] - np.einsum('ij,ij->i', ah[:, 1:], W[:, 1:])
    vv = (ah[:, 0:1] * W[:, 1:] + W[:, 0:1] * ah[:, 1:]
          + np.cross(ah[:, 1:], W[:, 1:]))
    v = np.concatenate([v0[:, None], vv], axis=1)
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    U[mu][idx] = quat_to_su2(v)


def sweep(U, L, beta, rng):
    ps = parity_sets(L)
    for c in (0, 1):
        for mu in range(4):
            heatbath_link(U, mu, ps[c], beta, rng)


def plaquettes(U):
    P = np.empty((6,) + U.shape[1:5] + (2, 2), dtype=complex)
    k = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            P[k] = U[mu] @ np.roll(U[nu], -1, axis=mu) @ dag(np.roll(U[mu], -1, axis=nu)) @ dag(U[nu])
            k += 1
    return P


def q_planes(P):
    """q = 1 - (1/2) Re tr U_plaq,  shape (6, L,L,L,L)"""
    return 1.0 - 0.5 * P.real.trace(axis1=-2, axis2=-1)


# ---------------------------------------------------------------- validation V4
def check_staple2(L=2, seed=1):
    """brute force: sum over the 6 plaquettes containing the link, all with U_mu(x) first."""
    rng = np.random.default_rng(seed)
    U = random_su2(rng, (4, L, L, L, L))
    ps = parity_sets(L)
    err = 0.0
    for mu in range(4):
        for c in (0, 1):
            idx = ps[c]
            S = staple_mu(U, mu, idx)
            lhs = (U[mu][idx] @ S).real.trace(axis1=-2, axis2=-1)
            rhs = np.zeros(len(idx[0]))
            for nu in range(4):
                if nu == mu:
                    continue
                # plaquette at x: U_mu(x), U_nu(x+mu), U_mu(x+nu)^dag, U_nu(x)^dag
                P1 = U[mu] @ np.roll(U[nu], -1, axis=mu) @ dag(np.roll(U[mu], -1, axis=nu)) @ dag(U[nu])
                rhs += P1[idx].real.trace(axis1=-2, axis2=-1)
                # plaquette at x-nu: contains U_mu(x-nu) and U_mu(x)^dag -> use reality
                xm = tuple(np.roll(idx[a], 1, axis=nu) for a in range(4))
                P2 = U[mu] @ np.roll(U[nu], -1, axis=mu) @ dag(np.roll(U[mu], -1, axis=nu)) @ dag(U[nu])
                rhs += P2[xm].real.trace(axis1=-2, axis2=-1)
            err = max(err, np.abs(lhs - rhs).max())
    return err


# ---------------------------------------------------------------- Metropolis (independent check)
def metropolis(U, L, beta, rng, nsweep, eps=0.5):
    ps = parity_sets(L)
    trP = []
    for it in range(nsweep):
        for c in (0, 1):
            for mu in range(4):
                idx = ps[c]
                S = staple_mu(U, mu, idx)
                old = U[mu][idx]
                # random small SU(2): exp(i eps * n.sigma)
                n = rng.standard_normal((len(idx[0]), 3))
                n /= np.linalg.norm(n, axis=1, keepdims=True)
                th = eps * rng.random(len(idx[0]))
                V = (np.cos(th)[:, None, None] * ID2
                     + 1j * np.sin(th)[:, None, None]
                     * (n[:, 0, None, None] * S1 + n[:, 1, None, None] * S2 + n[:, 2, None, None] * S3))
                new = V @ old
                dS = -0.5 * beta * ((new @ S).real.trace(axis1=-2, axis2=-1)
                                    - (old @ S).real.trace(axis1=-2, axis2=-1))
                acc = rng.random(len(idx[0])) < np.exp(np.minimum(dS, 0.0))
                U[mu][idx] = np.where(acc[:, None, None], new, old)
        P = plaquettes(U)
        trP.append(0.5 * P.real.trace(axis1=-2, axis2=-1).mean())
    return np.array(trP)


# ---------------------------------------------------------------- measurement
def measure_run(L, beta, ntherm, nkeep, nstreams, seed0=1000, gap=1):
    ps = parity_sets(L)
    V = L ** 4
    out = []
    t0 = time.time()
    for s in range(nstreams):
        rng = np.random.default_rng(seed0 + s)
        U = random_su2(rng, (4, L, L, L, L))
        for _ in range(ntherm):
            sweep(U, L, beta, rng)
        # accumulators
        n = 0
        sP = 0.0
        sq1 = 0.0
        sq1_sq = 0.0
        M1 = np.zeros((V, V))       # <q1(x) q1(y)>
        M6 = np.zeros((V, V))
        sum1 = np.zeros(V)
        sum6 = np.zeros(V)
        series = []
        for it in range(nkeep):
            sweep(U, L, beta, rng)
            if it % gap:
                continue
            q = q_planes(plaquettes(U))
            q1 = q[0].reshape(-1).real
            q6 = q.sum(axis=0).reshape(-1).real
            M1 += np.outer(q1, q1)
            M6 += np.outer(q6, q6)
            sum1 += q1
            sum6 += q6
            series.append(q1.mean())
            n += 1
            sP += (1.0 - q.mean())   # mean over 24 plaquettes of (1-q) = P
        M1 /= n; M6 /= n; sum1 /= n; sum6 /= n
        out.append(dict(n=n, M1=M1, M6=M6, mean1=sum1, mean6=sum6,
                        Pmean=sP / n, series=np.array(series)))
    return out, time.time() - t0


def var_block(M, mean, mask_idx):
    """Var of block average over site set given by mask_idx (list of site indices)."""
    nb = len(mask_idx)
    S = M[np.ix_(mask_idx, mask_idx)].sum()
    m = mean[mask_idx].sum()
    return S / nb ** 2 - (m / nb) ** 2


def analyse(streams, L):
    V = L ** 4
    site = np.arange(V)
    res = {}
    # per-stream estimates
    keys = ['P', 'varq1', 'varq6', 'varqbar1', 'varqbar6', 'chi1', 'chi6']
    per = {k: [] for k in keys}
    for st in streams:
        M1, M6 = st['M1'], st['M6']
        m1, m6 = st['mean1'], st['mean6']
        per['P'].append(st['Pmean'])
        per['varq1'].append(np.diag(M1 - np.outer(m1, m1)).mean())
        per['varq6'].append(np.diag(M6 - np.outer(m6, m6)).mean())
        v1 = var_block(M1, m1, site)
        v6 = var_block(M6, m6, site)
        per['varqbar1'].append(v1)
        per['varqbar6'].append(v6)
        per['chi1'].append(V * v1)
        per['chi6'].append(V * v6)
    ns = len(streams)
    for k in keys:
        a = np.array(per[k])
        mean = a.mean()
        se = a.std(ddof=1) / np.sqrt(ns) if ns > 1 else 0.0
        res[k] = (float(mean), float(se))
    # correlator shells for q1
    M1 = np.mean([s['M1'] for s in streams], axis=0)
    m1 = np.mean([s['mean1'] for s in streams], axis=0)
    C = M1 - np.outer(m1, m1)
    coords = np.array(np.unravel_index(np.arange(V), (L,) * 4)).T
    shells = {}
    for i in range(V):
        for j in range(i + 1, V):
            d = coords[i] - coords[j]
            d = np.where(d > L // 2, d - L, d)
            d = np.where(d < -L // 2, d + L, d)
            d = np.where(d == L // 2, d - L, d)   # L=2: +1 -> -1
            r = float(np.linalg.norm(d))
            shells.setdefault(r, []).append(C[i, j])
    shell_out = {f"r={r}": (float(np.mean(v)), int(len(v))) for r, v in sorted(shells.items())}
    # cross-check C(x,y) translation invariance
    dev = 0.0
    for i in range(V):
        for j in range(V):
            d = coords[i] - coords[j]
            d = (d % L)
            r = float(np.linalg.norm(np.where(d > L // 2, d - L, d)))
            dev = max(dev, abs(C[i, j] - shell_out[f"r={r}"][0]))
    res['_shells'] = shell_out
    res['_C_dev'] = float(dev)
    res['_tau'] = autocorr_tau(np.array([s['series'] for s in streams]))
    return res


def autocorr_tau(series_list, maxlag=400):
    taus = []
    for s in series_list:
        s = s - s.mean()
        var = np.mean(s * s)
        if var <= 0:
            continue
        ac = np.array([np.mean(s[:-k] * s[k:]) / var for k in range(1, min(maxlag, len(s) // 2))])
        # Sokal windowing
        tau = 1.0
        for k, v in enumerate(ac, start=1):
            tau += 2 * v
            if k >= 6 * max(tau, 1.0):
                break
        taus.append(tau)
    return float(np.max(taus)) if taus else float('nan')


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
    log = []

    def P(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        log.append(s)

    P("=" * 78)
    P("HZW-L2 : varpi = a^4 b^4 Z_a at L=2 SU(2) beta=2.5   (a=1 lattice units)")
    P("=" * 78)

    P("")
    P("[V4] staple identity Re tr(U S) = sum_6plaq Re tr(P)")
    P("     max|err| = %.3e" % check_staple2())

    P("")
    P("[V1] beta=0 anchor: Var(q) should be 1/4 exactly")
    rng = np.random.default_rng(7)
    U = random_su2(rng, (4, 2, 2, 2, 2))
    qs = []
    for _ in range(4000):
        sweep(U, 2, 0.0, rng)
        qs.append(q_planes(plaquettes(U))[0].reshape(-1))
    qs = np.concatenate(qs)
    P("     <q> = %.5f   Var(q) = %.5f   (exact 0.25)" % (qs.mean(), qs.var()))

    P("")
    P("[V2] L=4 beta=2.5 <P> (repo: 0.6526)")
    rng = np.random.default_rng(11)
    U = random_su2(rng, (4, 4, 4, 4, 4))
    for _ in range(400):
        sweep(U, 4, 2.5, rng)
    trP = []
    for _ in range(2000):
        sweep(U, 4, 2.5, rng)
        P4 = plaquettes(U)
        trP.append(0.5 * P4.real.trace(axis1=-2, axis2=-1).mean())
    P("     <P>_L4 = %.4f +- %.4f  (n=2000)" % (np.mean(trP), np.std(trP) / np.sqrt(len(trP))))

    P("")
    P("[V3] L=2 beta=2.5 : heat bath vs Metropolis")
    rng = np.random.default_rng(21)
    U = random_su2(rng, (4, 2, 2, 2, 2))
    m = metropolis(U, 2, 2.5, rng, 4000)
    P("     Metropolis <P> = %.4f +- %.4f" % (m[500:].mean(), m[500:].std() / np.sqrt(len(m[500:]))))

    # production
    L = 2
    beta = 2.5
    ntherm, nkeep, nstreams = 1200, 10000, 8
    t0 = time.time()
    streams, dt = measure_run(L, beta, ntherm, nkeep, nstreams, seed0=5000)
    res = analyse(streams, L)
    ncfg = sum(s['n'] for s in streams)
    P("")
    P("[1] production  L=2 beta=2.5  streams=%d  sweeps/stream=%d (therm %d)  configs=%d  wall=%.1fs"
      % (nstreams, nkeep, ntherm, ncfg, time.time() - t0))
    P("     tau_int(qbar) (max over streams) = %.2f ; N_eff = %d" % (res['_tau'], ncfg / max(res['_tau'], 1e-9)))
    P("     <P> = %.5f +- %.5f      (P = (1/2)Re tr U_plaq, single plane)" % res['P'])
    P("     Var(q) single plane  = %.5f +- %.5f" % res['varq1'])
    P("     Var(q6) six-plane sum= %.5f +- %.5f" % res['varq6'])
    P("     Var(qbar) block=L^4  = %.6f +- %.6f" % res['varqbar1'])
    P("     chi = V*Var(qbar)    = %.5f +- %.5f   (= varpi at b=L=2a, N_B=16)" % res['chi1'])
    P("     varpi(b=2a) = 16*Var(qbar) = %.5f +- %.5f" % res['chi1'])
    P("     varpi(b=a)  = Var(q)       = %.5f +- %.5f" % res['varq1'])
    P("     ratio varpi(2a)/varpi(a)   = %.3f" % (res['chi1'][0] / res['varq1'][0]))
    P("     six-plane: varpi6(b=2a) = %.5f +- %.5f ; varpi6(b=a) = %.5f +- %.5f"
      % (res['chi6'][0], res['chi6'][1], res['varq6'][0], res['varq6'][1]))
    P("     kinematic ceiling (R1) 4*N_B = 64  -> varpi/ceiling = %.2e" % (res['chi1'][0] / 64.0))
    P("")
    P("     C_q1 shells (|z|, mean, multiplicity):")
    for k, v in res['_shells'].items():
        P("        %-6s %+.5f   n=%d" % (k, v[0], v[1]))
    P("     translation-invariance max dev of C(x,y) = %.2e" % res['_C_dev'])

    js = dict(L=L, beta=beta, nstreams=nstreams, nkeep=nkeep, ncfg=ncfg,
              P=res['P'], varq1=res['varq1'], varq6=res['varq6'],
              varqbar1=res['varqbar1'], chi1=res['chi1'], chi6=res['chi6'],
              shells={k: v[0] for k, v in res['_shells'].items()},
              tau=res['_tau'], staple_err=float(check_staple2()))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hzw_l2_za.json'), 'w') as f:
        json.dump(js, f, indent=1)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hzw_l2_za.log'), 'w') as f:
        f.write("\n".join(log) + "\n")
    print("\n# log -> hzw_l2_za.log")
