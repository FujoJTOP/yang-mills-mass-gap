"""
W4 -- blocked / collective positivity probe.

Question: the 04-savvidy number (per-link Hessian ~49% negative at beta=2.5)
is a SINGLE-LINK, STAPLE-FROZEN, LEFT-FRAME curvature.  The object the SLC
hypothesis actually needs is the GAUGE-INVARIANT HESSIAN OF THE ACTION ON THE
COLLECTIVE / BLOCKED (low-momentum) TANGENT SPACE.

We build:
  (0) per-link staple-frozen Hessian (reproduce the 04 number)          -- baseline
  (1) the intrinsic manifold Hessian H = d^2 S on the full link tangent
      (left-geodesic second derivative; chart-independent bilinear form)
  (2) project out gauge orbits  -> transverse (Faddeev-Popov-slice) H_t
  (3) low-pass "blocked" subspace (k_mu <= pi/2, block factor 2) -> H_b
Then we report negative-eigenvalue fraction and minimum for each, at beta=2.5,
microscopic, and after APE smoothing (collective smoothing = flow/block proxy).

Conventions: SU(2) Wilson action S = beta * sum_p (1 - 1/2 Re tr U_p).
t^a = sigma^a/2, tr(t^a t^b) = 1/2 delta_ab.  Left variation dU = i eps t^a U.
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from su3 import su2_heatbath, project_su3, ID2, SIG  # noqa: E402

T2 = SIG / 2.0                      # (3,2,2) Hermitian generators


def adj(A):
    return np.conj(np.swapaxes(A, -1, -2))


def expm_t(a, h):
    """exp(i h t^a) for SU(2)."""
    return np.cos(h / 2.0) * ID2 + 1j * np.sin(h / 2.0) * SIG[a]


# ----------------------------------------------------------------------
# SU(2) lattice machinery
# ----------------------------------------------------------------------
def random_su2(rng, shape=()):
    z = (rng.standard_normal(shape + (2, 2)) +
         1j * rng.standard_normal(shape + (2, 2))) / np.sqrt(2.0)
    q, r = np.linalg.qr(z)
    d = np.diagonal(r, axis1=-2, axis2=-1)
    q = q * (d / np.abs(d))[..., None, :]
    det = np.linalg.det(q)
    q[..., 0, :] *= np.conj(det)[..., None]
    return q


def staple_sum2(U, mu):
    L = U.shape[1]
    Um = U[mu]
    S = np.zeros((L, L, L, L, 2, 2), dtype=np.complex128)
    for nu in range(4):
        if nu == mu:
            continue
        Un = U[nu]
        S1 = np.roll(Un, -1, axis=mu) @ adj(np.roll(Um, -1, axis=nu)) @ adj(Un)
        D = np.roll(np.roll(Un, 1, axis=nu), -1, axis=mu)
        E = np.roll(Um, 1, axis=nu)
        F = np.roll(Un, 1, axis=nu)
        S += adj(D) @ adj(E) @ F
        S += S1
    return S


def action2(U, beta):
    L = U.shape[1]
    tot = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Up = (U[mu] @ np.roll(U[nu], -1, axis=mu)
                  @ adj(np.roll(U[mu], -1, axis=nu)) @ adj(U[nu]))
            tot += np.real(np.trace(Up, axis1=-2, axis2=-1)).sum() / 2.0
    return beta * (6.0 * L ** 4 - tot)


def grad2(U, beta):
    """Left-frame gradient G[mu,x,a] = (beta/2) Im tr(t^a U S)."""
    G = np.empty((4,) + U.shape[1:5] + (3,))
    for mu in range(4):
        X = U[mu] @ staple_sum2(U, mu)
        G[mu] = (beta / 2.0) * np.einsum('aij,...ji->...a', T2, X).imag
    return G


def perlink_min2(U, beta):
    """Staple-frozen per-link Hessian min eigenvalue (the 04-style number)."""
    L = U.shape[1]
    out = np.empty((4, L, L, L, L))
    for mu in range(4):
        X = U[mu] @ staple_sum2(U, mu)
        H = np.einsum('aij,bjk,...ki->...ab', T2, T2, X).real
        H = 0.5 * (H + np.swapaxes(H, -1, -2))
        out[mu] = (beta / 2.0) * np.linalg.eigvalsh(H)[..., 0]
    return out


# ----------------------------------------------------------------------
# intrinsic Hessian (left-geodesic second derivative)
# ----------------------------------------------------------------------
def hessian2(U, beta, h=1e-5):
    L = U.shape[1]
    ns = L ** 4
    dim = 4 * ns * 3
    H = np.empty((dim, dim))
    E = {1: [expm_t(a, h) for a in range(3)],
         -1: [expm_t(a, -h) for a in range(3)]}
    base = U.copy()
    col = 0
    for mu in range(4):
        for s in range(ns):
            for a in range(3):
                Up = base.copy()
                Up.reshape(4, ns, 2, 2)[mu, s] = E[1][a] @ base.reshape(4, ns, 2, 2)[mu, s]
                Um = base.copy()
                Um.reshape(4, ns, 2, 2)[mu, s] = E[-1][a] @ base.reshape(4, ns, 2, 2)[mu, s]
                gp = grad2(Up, beta).reshape(-1)
                gm = grad2(Um, beta).reshape(-1)
                H[:, col] = (gp - gm) / (2.0 * h)
                col += 1
    return 0.5 * (H + H.T)


# ----------------------------------------------------------------------
# gauge orbit basis and transverse / low-pass projectors
# ----------------------------------------------------------------------
def _site_index_tables(L):
    ns = L ** 4
    idx4 = np.arange(ns).reshape(L, L, L, L)
    coords = np.stack(np.meshgrid(*[np.arange(L)] * 4, indexing='ij'),
                      axis=-1).reshape(ns, 4)
    return ns, idx4, coords


def comp2(X):
    """components v_c of Hermitian X in t^c basis: v_c = 2 Re tr(t^c X)."""
    return 2.0 * np.real(np.einsum('cij,ji->c', T2, X))


def gauge_basis2(U):
    """Orthonormal basis Q (dim x 3*ns) of gauge-orbit tangent (left frame)."""
    L = U.shape[1]
    ns, idx4, coords = _site_index_tables(L)
    dim = 4 * ns * 3

    def fidx(mu, s, a):
        return (mu * ns + s) * 3 + a

    Uf = U.reshape(4, ns, 2, 2)
    Q = np.zeros((dim, ns * 3))
    col = 0
    for s in range(ns):
        for a in range(3):
            vec = np.zeros(dim)
            for mu in range(4):
                # link (mu, s): + t^a
                vec[fidx(mu, s, a)] += 1.0
                # link (mu, y) with y = s - mu_hat : - U t^a U^dag
                cy = coords[s].copy()
                cy[mu] = (cy[mu] - 1) % L
                y = idx4[cy[0], cy[1], cy[2], cy[3]]
                V = -Uf[mu, y] @ T2[a] @ adj(Uf[mu, y])
                vin = comp2(V)
                for c in range(3):
                    vec[fidx(mu, y, c)] += vin[c]
            Q[:, col] = vec
            col += 1
    # orthonormalize
    Q, _ = np.linalg.qr(Q)
    return Q


def transverse_basis(Q):
    """Orthonormal basis Z (dim x (dim-ng)) of the gauge-orthogonal complement."""
    dim = Q.shape[0]
    ng = Q.shape[1]
    _, S, Vh = np.linalg.svd(Q.T)
    rank = int(np.sum(S > 1e-8))
    return Vh[rank:].T


def lowpass_basis(U, L, kfrac=0.5):
    """Complex Fourier low-pass basis (dim x n) with |k_mu| <= kfrac*pi
    (block factor 1/kfrac).  Includes both signs via full k grid."""
    ns = L ** 4
    dim = 4 * ns * 3
    k1dg = 2.0 * np.pi * np.arange(L) / L
    kmax = kfrac * np.pi + 1e-9
    keep = [kk for kk in k1dg if kk <= kmax]
    coords = np.stack(np.meshgrid(*[np.arange(L)] * 4, indexing='ij'),
                      axis=-1).reshape(ns, 4)
    cols = []
    from itertools import product
    for kvec in product(keep, repeat=4):
        phase = np.exp(1j * (coords @ np.array(kvec))) / np.sqrt(ns)
        for mu in range(4):
            for a in range(3):
                v = np.zeros(dim, dtype=np.complex128)
                v[(mu * ns + np.arange(ns)) * 3 + a] = phase
                cols.append(v)
    return np.array(cols).T


def orthonormalize_against(W, Q):
    """Remove span(Q) from columns of W, then orthonormalize."""
    Wp = W - Q @ (Q.conj().T @ W)
    U_, S, _ = np.linalg.svd(Wp, full_matrices=False)
    r = int(np.sum(S > 1e-8))
    return U_[:, :r]


# ----------------------------------------------------------------------
# measurement driver
# ----------------------------------------------------------------------
def sample2(L, beta, rng, ntherm, nkeep, gapevery=5):
    U = random_su2(rng, (4, L, L, L, L))
    ns = L ** 4
    idx = np.arange(ns).reshape(L, L, L, L)
    x0, x1, x2, x3 = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
    par = (x0 + x1 + x2 + x3) % 2
    idxp = [np.nonzero(par.reshape(-1) == c)[0] for c in (0, 1)]
    kappa = beta / 2.0

    def sweep(U):
        Uf = U.reshape(4, ns, 2, 2)
        for c in (0, 1):
            ids = idxp[c]
            for mu in range(4):
                S = staple_sum2(U, mu).reshape(ns, 2, 2)[ids]
                A = Uf[mu, ids]
                M = A @ S
                v = su2_heatbath(kappa, M, 0, 1, rng)
                Uf[mu, ids] = v @ A

    for _ in range(ntherm):
        sweep(U)
    out = []
    for k in range(nkeep):
        for _ in range(gapevery):
            sweep(U)
        out.append(U.copy())
    return out


def ape_all2(U, alpha=0.5, nsteps=1):
    """Isotropic gauge-covariant APE smoothing (all 4 directions) = block proxy."""
    U = U.copy()
    for _ in range(nsteps):
        old = U.copy()
        for mu in range(4):
            S = np.zeros_like(U[mu])
            for nu in range(4):
                if nu == mu:
                    continue
                # forward staple
                S += old[nu] @ np.roll(old[mu], -1, axis=nu) @ adj(np.roll(old[nu], -1, axis=mu))
                # backward staple
                S += (adj(np.roll(old[nu], 1, axis=nu)) @ np.roll(old[mu], 1, axis=nu)
                      @ np.roll(np.roll(old[nu], 1, axis=nu), -1, axis=mu))
            U[mu] = project_su3((1.0 - alpha) * old[mu] + alpha / 6.0 * S)
    return U


def decimate2(U):
    """One decimation blocking L -> L/2: U'_mu(2x) = proj( U_mu(2x) U_mu(2x+mu) )."""
    L = U.shape[1]
    Lp = L // 2
    idx = np.arange(0, L, 2)
    Ub = np.empty((4, Lp, Lp, Lp, Lp, 2, 2), dtype=np.complex128)
    for mu in range(4):
        A = U[mu]
        B = np.roll(A, -1, axis=mu)
        W = A @ B
        Ub[mu] = project_su3(W[np.ix_(idx, idx, idx, idx)])
    return Ub


def analyse(U, beta, tag, do_lowpass=True, do_ape=False):
    L = U.shape[1]
    ns = L ** 4
    dim = 4 * ns * 3
    out = {}

    # baseline per-link
    pl = perlink_min2(U, beta)
    out['perlink_fracneg'] = float(np.mean(pl < 0))
    out['perlink_min_over_beta'] = float(pl.min() / beta)
    out['perlink_mean_over_beta'] = float(pl.mean() / beta)

    H = hessian2(U, beta)
    Q = gauge_basis2(U)
    ng = Q.shape[1]
    # validation: gradient orthogonal to gauge orbit
    g = grad2(U, beta).reshape(-1)
    out['grad_gauge_orth'] = float(np.max(np.abs(Q.T @ g)))

    Z = transverse_basis(Q)
    Ht = Z.T @ H @ Z
    Ht = 0.5 * (Ht + Ht.T)
    ev = np.linalg.eigvalsh(Ht)
    out['trans_dim'] = int(ev.size)
    out['trans_fracneg'] = float(np.mean(ev < 0))
    out['trans_min_over_beta'] = float(ev.min() / beta)
    out['trans_mean_over_beta'] = float(ev.mean() / beta)
    out['trans_neg_frac_of_abs'] = float(
        -np.sum(ev[ev < 0]) / max(np.sum(np.abs(ev)), 1e-30))

    if do_lowpass:
        W = lowpass_basis(U, L, kfrac=0.5)
        Wt = orthonormalize_against(W, Q)
        Hb = Wt.conj().T @ H @ Wt
        Hb = 0.5 * (Hb + Hb.conj().T)
        evb = np.linalg.eigvalsh(Hb)
        out['block_dim'] = int(evb.size)
        out['block_fracneg'] = float(np.mean(evb < 0))
        out['block_min_over_beta'] = float(evb.min() / beta)
        out['block_mean_over_beta'] = float(evb.mean() / beta)

    if do_ape:
        # low-momentum weight of each full-transverse eigenvector, to say whether
        # negative modes are IR or UV.
        w, V = np.linalg.eigh(Ht)
        # low-pass projector on transverse space
        W = lowpass_basis(U, L, kfrac=0.5)
        Wt = orthonormalize_against(W, Q)
        # weight of mode = || Wt^dag (Z @ v) ||^2  ; Z@v lives in full space
        # precompute G = Wt^dag Z  (low x trans)
        G = Wt.conj().T @ Z
        lowfrac = np.sum(np.abs(G @ V) ** 2, axis=0)  # per eigenvector
        neg = w < 0
        denom = np.sum(-w[neg]) if np.any(neg) else 1.0
        out['neg_density_low_frac'] = float(
            np.sum((-w[neg]) * lowfrac[neg]) / max(denom, 1e-30))
        out['neg_density_high_frac'] = float(
            np.sum((-w[neg]) * (1.0 - lowfrac[neg])) / max(denom, 1e-30))

    out['tag'] = tag
    return out


def run_case(L, beta, ntherm, nkeep, gapevery, seed, do_ape, do_decimate):
    rng = np.random.default_rng(seed)
    print(f"# sampling SU(2) L={L} beta={beta} ntherm={ntherm} nkeep={nkeep}", flush=True)
    cfgs = sample2(L, beta, rng, ntherm, nkeep, gapevery)
    rows = []
    for i, U in enumerate(cfgs):
        r = analyse(U, beta, f"cfg{i}", do_lowpass=True, do_ape=True)
        rows.append(r)
        print(f"  cfg{i}: perlink_neg={r['perlink_fracneg']:.3f} "
              f"trans_neg={r['trans_fracneg']:.3f} "
              f"trans_min/b={r['trans_min_over_beta']:+.4f} "
              f"block_neg={r['block_fracneg']:.3f} "
              f"block_min/b={r['block_min_over_beta']:+.4f} "
              f"neglow={r['neg_density_low_frac']:.3f}", flush=True)
        if do_ape:
            Us = ape_all2(U, alpha=0.5, nsteps=1)
            ra = analyse(Us, beta, f"cfg{i}_ape1", do_lowpass=True, do_ape=True)
            print(f"    APE1: perlink_neg={ra['perlink_fracneg']:.3f} "
                  f"trans_neg={ra['trans_fracneg']:.3f} "
                  f"trans_min/b={ra['trans_min_over_beta']:+.4f} "
                  f"block_neg={ra['block_fracneg']:.3f} "
                  f"block_min/b={ra['block_min_over_beta']:+.4f}", flush=True)
            rows.append(ra)
        if do_decimate and L % 2 == 0:
            Ud = decimate2(U)
            rd = analyse(Ud, beta, f"cfg{i}_dec", do_lowpass=False, do_ape=False)
            print(f"    DEC L={Ud.shape[1]}: perlink_neg={rd['perlink_fracneg']:.3f} "
                  f"trans_neg={rd['trans_fracneg']:.3f} "
                  f"trans_min/b={rd['trans_min_over_beta']:+.4f}", flush=True)
            rows.append(rd)
    return rows


def validate(beta=2.5, L=3, seed=0):
    """Sanity: analytic gradient vs FD; gauge orthogonality."""
    rng = np.random.default_rng(seed)
    U = random_su2(rng, (4, L, L, L, L))
    b = beta
    g = grad2(U, b)
    h = 1e-6
    err = 0.0
    for mu in range(4):
        for s in range(L ** 4):
            for a in range(3):
                Up = U.copy(); Up.reshape(4, L ** 4, 2, 2)[mu, s] = \
                    expm_t(a, h) @ U.reshape(4, L ** 4, 2, 2)[mu, s]
                Um = U.copy(); Um.reshape(4, L ** 4, 2, 2)[mu, s] = \
                    expm_t(a, -h) @ U.reshape(4, L ** 4, 2, 2)[mu, s]
                fd = (action2(Up, b) - action2(Um, b)) / (2 * h)
                err = max(err, abs(fd - g.reshape(-1)[(mu * L ** 4 + s) * 3 + a]))
    print(f"[validate] max |analytic_grad - FD_grad| = {err:.3e}")
    Q = gauge_basis2(U)
    go = np.max(np.abs(Q.T @ g.reshape(-1)))
    print(f"[validate] max |gauge . gradient| = {go:.3e}")
    # stability at weak coupling: transverse Hessian must be positive
    return err, go


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', default='validate')
    ap.add_argument('--L', type=int, default=4)
    ap.add_argument('--beta', type=float, default=2.5)
    ap.add_argument('--ntherm', type=int, default=300)
    ap.add_argument('--nkeep', type=int, default=4)
    ap.add_argument('--gapevery', type=int, default=5)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--ape', action='store_true')
    ap.add_argument('--decimate', action='store_true')
    args = ap.parse_args()
    if args.mode == 'validate':
        validate(beta=args.beta, L=args.L, seed=args.seed)
    else:
        run_case(args.L, args.beta, args.ntherm, args.nkeep, args.gapevery,
                 args.seed, args.ape, args.decimate)
