"""NA-3-T: Dobrushin / SDPI coefficients for the SU(3) single-link heat bath.

Only observables: the dynamics (su3.Lattice / su3_numba) is untouched.

Definitions (NA-3-T-spec.md §1):
  P_e(.|U_{≠e}) ∝ exp( (β/3) Re tr(U_e S_e) ) dU,   S_e = staple sum
  I(e,e') = sup_{U,U' differ only at e'} ||P_e(.|U) - P_e(.|U')||_TV
  η = max_e Σ_{e' ~ e} I(e,e')

Estimators here
  * rigorous upper bound: TV(p,q) ≤ (κ/2) ||ΔS||_*  (κ=β/3, nuclear norm)
  * numerical TV: Bennett acceptance-ratio (BAR) estimate of the normalisation
    ratio, then TV = 1 - ∫min(p,q).  Uses no cluster expansion.
  * No use of any defect/cluster weight.
"""
import numpy as np

import su3 as S

KAPPA = None  # set per call


# ----------------------------------------------------------------------
# conditional sampler: p(U) ∝ exp(κ Re tr(U Smat))
# ----------------------------------------------------------------------
def sample_cond(Smat, kappa, nsamp, nsweep, rng):
    U = S.random_su3(rng, (nsamp,))
    for _ in range(nsweep):
        M = U @ Smat
        for (i, j) in ((0, 1), (0, 2), (1, 2)):
            v = S.su2_heatbath(kappa, M, i, j, rng)
            Mi = M[:, i, :].copy(); Mj = M[:, j, :].copy()
            M[:, i, :] = v[:, 0, 0, None] * Mi + v[:, 0, 1, None] * Mj
            M[:, j, :] = v[:, 1, 0, None] * Mi + v[:, 1, 1, None] * Mj
            Ui = U[:, i, :].copy(); Uj = U[:, j, :].copy()
            U[:, i, :] = v[:, 0, 0, None] * Ui + v[:, 0, 1, None] * Uj
            U[:, j, :] = v[:, 1, 0, None] * Ui + v[:, 1, 1, None] * Uj
    return U


def _ftr(U, Smat):
    return np.real(np.einsum('nij,ji->n', U, Smat))


def tv_bar(Sp, Sq, kappa, nsamp=120, nsweep=6, rng=None):
    """BAR estimate of TV(p,q) for p∝e^{κRe tr(U Sp)}, q∝e^{κRe tr(U Sq)}."""
    rng = rng or np.random.default_rng(0)
    Up = sample_cond(Sp, kappa, nsamp, nsweep, rng)
    Uq = sample_cond(Sq, kappa, nsamp, nsweep, rng)
    a_p = kappa * _ftr(Up, Sp)
    b_p = kappa * _ftr(Up, Sq)
    a_q = kappa * _ftr(Uq, Sp)
    b_q = kappa * _ftr(Uq, Sq)
    # BAR: find L = log(Zp/Zq) solving
    #   mean_p 1/(1+exp((b-a) + log(Zq/Zp))) = mean_q 1/(1+exp((a-b) + log(Zp/Zq)))
    # with log(Zq/Zp) = -L.
    def fun(L):
        f = lambda x: 1.0 / (1.0 + np.exp(np.clip(x, -700, 700)))
        return np.mean(f(b_p - a_p - L)) - np.mean(f(a_q - b_q + L))
    lo, hi = -60.0, 60.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if fun(lo) * fun(mid) <= 0:
            hi = mid
        else:
            lo = mid
    L = 0.5 * (lo + hi)
    # ∫min(p,q) from the p-samples: min(1, q/p) = min(1, exp(b-a+L))
    rat = np.exp(np.clip(b_p - a_p + L, -700, 700))
    ov = np.mean(np.minimum(1.0, rat))
    return float(np.clip(1.0 - ov, 0.0, 1.0))


# ----------------------------------------------------------------------
# staple and neighbour structure
# ----------------------------------------------------------------------
def staple_at(U, x, mu):
    L = U.shape[1]
    Ssum = np.zeros((3, 3), dtype=np.complex128)
    for nu in range(4):
        if nu == mu:
            continue
        Un = U[nu]
        Um = U[mu]
        def at(m, y):
            return m[tuple(np.mod(y, L))]
        e = np.eye(4, dtype=int)
        yy = np.array(x) + e[mu]
        zz = np.array(x) + e[nu]
        Ssum += at(Un, yy) @ S.adj(at(Um, zz)) @ S.adj(at(Un, np.array(x)))
        w = np.array(x) - e[nu]
        v = w + e[mu]
        Ssum += S.adj(at(Un, v)) @ S.adj(at(Um, w)) @ at(Un, w)
    return Ssum


def neighbor_links(x, mu, L=None):
    """The links whose change can modify S_e (18 of them)."""
    e = np.eye(4, dtype=int)
    out = []
    for nu in range(4):
        if nu == mu:
            continue
        out.append((tuple(np.array(x) + e[nu]), mu))      # U_mu(x+nu)  in S1
        out.append((tuple(np.array(x) - e[nu]), mu))      # U_mu(x-nu)  in S2
        out.append((tuple(x), nu))                        # U_nu(x)
        out.append((tuple(np.array(x) + e[mu]), nu))      # U_nu(x+mu)
        out.append((tuple(np.array(x) - e[nu]), nu))      # U_nu(x-nu)
        out.append((tuple(np.array(x) - e[nu] + e[mu]), nu))
    if L is not None:
        out = [(tuple(int(c) % L for c in xx), mm) for (xx, mm) in out]
    # deduplicate
    seen, uniq = set(), []
    for item in out:
        if item not in seen:
            seen.add(item)
            uniq.append(item)
    return uniq


def nuclear_norm(A):
    return float(np.sum(np.linalg.svd(A, compute_uv=False)))


def tv_witness_lb(Sp, Sq, kappa, nsamp=400, nsweep=12, nsigma=4.0, rng=None):
    """Rigorous-in-the-statistical-sense lower bound on TV(p,q).

    Witness f(U) = Re tr(U (Sq-Sp)) / ||Sq-Sp||_*  in [-1,1]  =>
        TV >= (1/2) |E_p f - E_q f|.
    Returns (lb, value, se) with lb = value - nsigma*se, clipped at 0.
    """
    rng = rng or np.random.default_rng(0)
    dS = Sq - Sp
    nd = nuclear_norm(dS)
    if nd < 1e-14:
        return 0.0, 0.0, 0.0
    Up = sample_cond(Sp, kappa, nsamp, nsweep, rng)
    Uq = sample_cond(Sq, kappa, nsamp, nsweep, rng)
    fp = _ftr(Up, dS) / nd
    fq = _ftr(Uq, dS) / nd
    mp, sp = fp.mean(), fp.std(ddof=1) / np.sqrt(nsamp)
    mq, sq = fq.mean(), fq.std(ddof=1) / np.sqrt(nsamp)
    val = 0.5 * abs(mp - mq)
    se = 0.5 * np.sqrt(sp * sp + sq * sq)
    return float(max(0.0, val - nsigma * se)), float(val), float(se)


# ----------------------------------------------------------------------
# exact i.i.d. sampling of the link conditional (rejection from Haar)
# ----------------------------------------------------------------------
def sample_cond_exact(Smat, kappa, n, rng, chunk=100000):
    """Exact i.i.d. draws from p(U) ∝ exp(kappa * Re tr(U Smat)).

    Envelope: max_U Re tr(U Smat) = ||Smat||_* (nuclear norm), so Haar
    proposals are accepted with probability exp(kappa (Re tr(US) - ||S||_*)).
    Exact by construction; acceptance rate = Z e^{-kappa ||S||_*}.
    """
    M = kappa * nuclear_norm(Smat)
    out, got = [], 0
    while got < n:
        U = S.random_su3(rng, (chunk,))
        acc = np.exp(kappa * _ftr(U, Smat) - M)
        keep = U[rng.random(chunk) < acc]
        out.append(keep)
        got += len(keep)
    return np.concatenate(out)[:n]


def tv_ks(Sp, Sq, kappa, n=10000, delta=1e-6, rng=None):
    """EXACT TV between the two matrix-Fisher conditionals, by the KS gap.

    p/q = (Zq/Zp) exp(kappa Re tr(U (Sp-Sq))), a monotone function of the
    single scalar  T = Re tr(U (Sp-Sq)) / ||Sp-Sq||_* .  Hence {p>q} is a
    threshold set {T > c} for some c, and
        TV = E_p[1{T>c}] - E_q[1{T>c}] = sup_c (F_p(c) - F_q(c))
    up to orientation.  No BAR, no coupling, no cluster weight.
    Error: two-sample DKW, distribution-free: |F_n - F| <= sqrt(ln(2/delta)/(2n)).
    Returns (lb, value, err) with lb = value - 2*err.
    """
    rng = rng or np.random.default_rng(0)
    dS = Sp - Sq
    nd = nuclear_norm(dS)
    if nd < 1e-14:
        return 0.0, 0.0, 0.0
    Up = sample_cond_exact(Sp, kappa, n, rng)
    Uq = sample_cond_exact(Sq, kappa, n, rng)
    Tp = np.sort(_ftr(Up, dS) / nd)
    Tq = np.sort(_ftr(Uq, dS) / nd)
    grid = np.concatenate([Tp, Tq])
    Fp = np.searchsorted(Tp, grid, side='right') / n
    Fq = np.searchsorted(Tq, grid, side='right') / n
    val = float(max((Fp - Fq).max(), (Fq - Fp).max()))
    err = float(np.sqrt(np.log(2.0 / delta) / (2.0 * n)))
    return float(max(0.0, val - 2.0 * err)), val, err
