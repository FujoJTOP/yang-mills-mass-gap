"""
RERUN core: exact marginal curvature G''(0) with the CORRECT decomposition

    G''(0) = E[P^T H P | a=0] - Var(P^T grad S | a=0) + W

in FLAT (exponential-chart) coordinates.  We work in the chart  U = exp(i X),
X = sum_a x_a t^a Hermitian; the collective mode is a fixed unit algebra vector
v, and  a = v . x  (= v^T log U).  Flat displacements are x -> x + s v  i.e.
X -> X + s V  with V = sum_a v_a t^a  (NOT the group geodesic exp(isV)U).

Everything is consistent in this chart:
  H_vv   = d^2/ds^2 S( exp(i(X+sV)) )|_0                 (action Hessian, flat)
  G_v    = d/ds   S( exp(i(X+sV)) )|_0                   (action force,   flat)
  Phi    = -log J = sum_links sum_{i<j} g(lam_i-lam_j)   (Haar Jacobian, flat)
  g(d)   = -2 log( sin(d/2)/(d/2) )                      (one pair)
Then EXACTLY (flat-coordinate constrained free energy):
  V_eff''(0) = E[ H_vv + Phi'' | a=0 ] - Var( G_v + Phi' | a=0 )
             = ( E[H_vv]-Var(G_v) ) + W
  W = E[Phi''] - Var(Phi') - 2 Cov(G_v, Phi')
The gold-standard  kappa_hist = curvature of -log P(a)  must equal RHS.
"""
import sys
import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/unconditional/w7')
import w4_blocked as W        # noqa: E402
import gpp_core as C          # noqa: E402   (su2_log_x, plane_mode)


# ----------------------------------------------------------------------
# Jacobian pair function  g(d) = -2 log(sin(d/2)/(d/2))   (d in radians)
# ----------------------------------------------------------------------
def g_pair(d):
    """vectorised; d = eigenvalue difference.  Handles d~0 and |d|>2pi robustly."""
    d = np.asarray(d, dtype=float)
    half = 0.5 * d
    out = np.empty_like(d)
    small = np.abs(d) < 1e-6
    # series  g = d^2/12 + d^4/1440 + ...
    ds = d[small]
    out[small] = ds * ds / 12.0 + ds ** 4 / 1440.0
    big = ~small
    db = d[big]
    s = np.sin(0.5 * db) / (0.5 * db)
    out[big] = -2.0 * np.log(np.abs(s))
    return out


def phi_from_eigs(lam):
    """lam: (...,N) eigenvalues of X. Returns Phi=-log J = sum_{i<j} g(lam_i-lam_j)."""
    N = lam.shape[-1]
    tot = np.zeros(lam.shape[:-1])
    for i in range(N):
        for j in range(i + 1, N):
            tot = tot + g_pair(lam[..., i] - lam[..., j])
    return tot


# ----------------------------------------------------------------------
# SU(2) / SU(3) X-matrix helpers
# ----------------------------------------------------------------------
def X_from_x_su2(x):
    """x:(...,3) -> Hermitian X = sum x_a t^a = 0.5 x.sigma, shape (...,2,2)."""
    return np.einsum('...a,aij->...ij', x, W.T2)


_TG = None


def X_from_x_su3(x):
    """x:(...,8) -> Hermitian X = sum x_a T^a (gell-mann/2), shape (...,3,3)."""
    global _TG
    if _TG is None:
        import su3
        _TG = su3.T_GEN
    return np.einsum('...a,aij->...ij', x, _TG)


def su3_log_herm(U):
    """Traceless Hermitian X (per link) with expi(X)=U exactly, handling the
    center twist: if sum(principal eigenangles)=+-2pi, shift the angle nearest 0
    by -+2pi so the sum is 0 and expi(X)=U (not the Z3-twisted U).  Also returns
    x_a = 2 tr(t^a X).  This is the branch-safe SU(3) log."""
    global _TG
    if _TG is None:
        import su3
        _TG = su3.T_GEN
    shape = U.shape[:-2]
    Uf = U.reshape(-1, 3, 3)
    w, V = np.linalg.eig(Uf)
    phi = np.angle(w)                                    # (M,3)
    s = 2 * np.pi * np.round(phi.sum(axis=1) / (2 * np.pi))
    need = np.abs(s) > 1e-9
    if np.any(need):
        k = np.argmin(np.abs(phi), axis=1)
        m = np.nonzero(need)[0]
        phi[m, k[m]] -= s[m]
    X = np.einsum('mij,mj,mkj->mik', V, phi, np.conj(V))
    x = 2.0 * np.real(np.einsum('aij,mji->ma', _TG, X))
    return X.reshape(shape + (3, 3)), x.reshape(shape + (8,))


def expi(X):
    """exp(i X) for Hermitian X (...,n,n)."""
    w, V = np.linalg.eigh(X)
    E = np.einsum('...ik,...k,...jk->...ij', V, np.exp(1j * w), np.conj(V))
    return E


def S_of_x(x, beta, Xfrom, action):
    U = expi(Xfrom(x))
    return action(U, beta)


def flat_derivs(x, Vdir, beta, Xfrom, action, h=1e-3):
    """flat d/ds and d^2/ds^2 of S(x + s Vdir) at s=0 via central differences."""
    Sp = S_of_x(x + h * Vdir, beta, Xfrom, action)
    S0 = S_of_x(x, beta, Xfrom, action)
    Sm = S_of_x(x - h * Vdir, beta, Xfrom, action)
    return (Sp - Sm) / (2 * h), (Sp - 2 * S0 + Sm) / (h * h)


def flat_phi_derivs(x, Vdir, Xfrom, h=1e-3):
    """flat d/ds and d^2/ds^2 of Phi = -log J  along x+s Vdir at s=0 (FD)."""
    def phi_of(xx):
        X = Xfrom(xx)
        lam = np.linalg.eigvalsh(X)          # (...,N)
        return float(np.sum(phi_from_eigs(lam)))   # sum over links
    Pp = phi_of(x + h * Vdir)
    P0 = phi_of(x)
    Pm = phi_of(x - h * Vdir)
    return (Pp - Pm) / (2 * h), (Pp - 2 * P0 + Pm) / (h * h)


def _gp(d):
    """g'(d) = -cot(d/2) + 2/d, with series."""
    d = np.asarray(d, float); out = np.empty_like(d)
    small = np.abs(d) < 1e-4; ds = d[small]
    out[small] = ds / 6.0 + ds ** 3 / 180.0
    b = ~small; db = d[b]
    out[b] = -1.0 / np.tan(0.5 * db) + 2.0 / db
    return out


def _gpp(d):
    """g''(d) = 0.5/sin^2(d/2) - 2/d^2, with series."""
    d = np.asarray(d, float); out = np.empty_like(d)
    small = np.abs(d) < 1e-4; ds = d[small]
    out[small] = 1.0 / 6.0 + ds ** 2 / 60.0
    b = ~small; db = d[b]
    out[b] = 0.5 / np.sin(0.5 * db) ** 2 - 2.0 / db ** 2
    return out


def phi_derivs_analytic(x, Vdir, Xfrom):
    """Analytic flat d/ds, d^2/ds^2 of Phi=-log J along x+s Vdir at s=0,
    via eigenvalue perturbation theory (no FD across branch cuts)."""
    X = Xfrom(x)
    Phi_ = Xfrom(Vdir)                        # displacement matrix V
    Xf = X.reshape(-1, X.shape[-2], X.shape[-1])
    Vf = Phi_.reshape(-1, X.shape[-2], X.shape[-1])
    w, U = np.linalg.eigh(Xf)                 # w (M,N), U (M,N,N)
    Uh = np.conj(np.swapaxes(U, -1, -2))
    Vm = Uh @ Vf @ U                          # eigenbasis (M,N,N)
    N = w.shape[-1]
    P1 = np.zeros(w.shape[0]); P2 = np.zeros(w.shape[0])
    d1 = np.real(np.diagonal(Vm, axis1=-2, axis2=-1))      # delta_i (M,N)
    # second-order eigenvalue shift  d2phi_i = 2 sum_{k!=i}|Vm_ki|^2/(phi_i-phi_k)
    d2 = np.zeros_like(d1)
    for i in range(N):
        for k in range(N):
            if k == i:
                continue
            dd = w[:, i] - w[:, k]
            d2[:, i] += 2.0 * np.abs(Vm[:, k, i]) ** 2 / np.where(np.abs(dd) < 1e-12, 1e-12, dd)
    for i in range(N):
        for j in range(i + 1, N):
            dij = w[:, i] - w[:, j]
            P1 += _gp(dij) * (d1[:, i] - d1[:, j])
            P2 += (_gpp(dij) * (d1[:, i] - d1[:, j]) ** 2
                   + _gp(dij) * (d2[:, i] - d2[:, j]))
    return float(P1.sum()), float(P2.sum())


# ----------------------------------------------------------------------
# curvature of -log P(a) (gold standard), same as W7's veff_curvature
# ----------------------------------------------------------------------
def kappa_kde(a, rng=None, nboot=200, bw_method=0.2, rfit=0.8):
    """V_eff''(0) = curvature of -log P(a) at the peak, via Gaussian KDE
    (quadratic fit of -log KDE over |a|<rfit*sd).  Bootstrap SE."""
    from scipy.stats import gaussian_kde
    a = np.asarray(a, float)
    if rng is None:
        rng = np.random.default_rng(0)

    def est(sample):
        sd = sample.std()
        kde = gaussian_kde(sample, bw_method=bw_method)
        xs = np.linspace(-rfit * sd, rfit * sd, 48)
        return 2.0 * np.polyfit(xs, -np.log(kde(xs)), 2)[0]

    k = est(a)
    boots = np.empty(nboot)
    for b in range(nboot):
        boots[b] = est(a[rng.integers(0, a.size, a.size)])
    return float(k), float(boots.std())


def kernel_cond(a, cols, bw=0.3):
    """Nadaraya-Watson conditional moments at a=0 with Gaussian kernel.
    cols: dict name->array(N,).  Returns dict of weighted E and pairwise Cov.
    Weighted variance/cov use the weighted (population) form."""
    a = np.asarray(a, float)
    sd = a.std()
    w = np.exp(-0.5 * (a / (bw * sd)) ** 2)
    ws = w.sum()
    out = {'n_eff': float(ws ** 2 / np.sum(w ** 2))}
    means = {}
    for k, x in cols.items():
        x = np.asarray(x, float)
        means[k] = float(np.sum(w * x) / ws)
        out['E_' + k] = means[k]
    for k, x in cols.items():
        for l, y in cols.items():
            if k <= l:
                xm = np.asarray(x, float) - means[k]
                ym = np.asarray(y, float) - means[l]
                out['C_' + k + '_' + l] = float(np.sum(w * xm * ym) / ws)
    out['V_own'] = {k: out['C_' + k + '_' + k] for k in cols}
    return out


def kappa_hist(a, nboot=300, rng=None, nave=2.5, nbins=40):
    if rng is None:
        rng = np.random.default_rng(0)

    def fit(sample):
        sd = sample.std()
        m = np.abs(sample) < nave * sd
        if m.sum() < nbins:
            m = np.abs(sample) < 4 * sd
        lo, hi = sample[m].min(), sample[m].max()
        edges = np.linspace(lo, hi, nbins + 1)
        cnt, _ = np.histogram(sample[m], bins=edges)
        cen = 0.5 * (edges[1:] + edges[:-1])
        good = cnt > 0
        cc = cen[good]
        yy = -np.log(cnt[good] / (sample[m].size * (edges[1] - edges[0])))
        A = np.vstack([np.ones_like(cc), cc ** 2, cc ** 4]).T
        coef, *_ = np.linalg.lstsq(A, yy, rcond=None)
        return 2.0 * coef[1], 4.0 * coef[2]

    kap, lam = fit(a)
    boots = np.empty(nboot)
    n = a.size
    for b in range(nboot):
        idx = rng.integers(0, n, n)
        boots[b] = fit(a[idx])[0]
    return kap, boots.std(), lam, float(1.0 / a.var())


if __name__ == '__main__':
    # ---- unit checks ----
    rng = np.random.default_rng(0)
    L = 3
    U = W.random_su2(rng, (4, L, L, L, L))
    x = C.su2_log_x(U)
    X = X_from_x_su2(x)
    U2 = expi(X)
    print("SU2  ||expi(X)-U||_max =", np.abs(U2 - U).max())
    xb = np.random.default_rng(1).standard_normal((5, 3))
    Xb = X_from_x_su2(xb)
    Ub = expi(Xb)
    xb2 = C.su2_log_x(Ub)
    print("SU2  log(expi(X))-x_max =", np.abs(xb2 - xb).max())
    # Jacobian: FD second derivative vs closed form w''(theta) etc.
    v = rng.standard_normal((4, L, L, L, L, 3))
    v /= np.linalg.norm(v)
    sig = 2.5
    gd, g2 = flat_derivs(x, v, sig, X_from_x_su2, W.action2)
    fp, fpp = flat_phi_derivs(x, v, X_from_x_su2)
    # closed-form Phi'' along v (per link) using w=g, theta=|x|
    th = np.linalg.norm(x, axis=-1)
    xhatv = np.sum(x * v, axis=-1)                # xhat . v_link times theta? careful
    # v_link . xhat = (x.v_link)/theta
    xv = np.sum(x * v, axis=-1)                   # x . v_link
    v2 = np.sum(v * v, axis=-1)
    # g(x)=g(|x|): Hessian directional = g''(th)*(xhat.v)^2 + (g'/th)*(|v|^2-(xhat.v)^2)
    def gpp(th):
        # g= -2log(sin(th/2)/(th/2)) ; g'' ; use series for small
        out = np.empty_like(th); sm = np.abs(th) < 1e-4
        out[sm] = 1.0/6.0 - 2.0*th[sm]**2/60.0
        t = th[~sm]
        out[~sm] = (0.5/np.sin(0.5*t)**2) - 2.0/t**2
        return out
    def gp_over_t(th):
        out = np.empty_like(th); sm = np.abs(th) < 1e-4
        out[sm] = 1.0/6.0 - th[sm]**2/60.0
        t = th[~sm]
        out[~sm] = (-1.0/np.tan(0.5*t) + 2.0/t)/t
        return out
    xhatv2 = np.where(th > 0, xv**2/np.maximum(th,1e-12)**2, 0.0)
    closed = np.sum(gpp(th)*xhatv2 + gp_over_t(th)*(v2 - xhatv2))
    print("SU2  Phi'' FD =%.6f  closed =%.6f" % (fpp, closed))
