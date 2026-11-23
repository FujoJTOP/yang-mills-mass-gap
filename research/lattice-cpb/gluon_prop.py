"""
Gluon propagator D(k) on SU(2) Wilson lattice in Landau gauge.

Wilson action S = beta * sum_p (1 - 1/2 Re tr U_p), beta = 4/g^2 (SU(2)), beta_c ~ 2.3.
Heat bath: w4_blocked.sample2 (kappa = beta/2, validated against Metropolis).

For each thermalized config:
 1. Landau gauge fix (Gauss-Seidel by parity colour, local optimal SU(2) update)
 2. gauge field  A^a_mu(x) = 2 Im tr( U_mu(x) t^a ),  t^a = sigma^a/2
 3. FFT -> At_mu^a(k);  G_munu(k) = (1/V) sum_a At_mu conj(At_nu)
 4. lattice transverse projector P^T_munu = delta - khat_mu khat_nu/khat^2,
    khat_mu = 2 sin(pi n_mu / L)
    D(k) = (1/3) Tr(P^T G)      (one transverse gluon, 3 colours averaged)
    DL(k)= (khat.khat/khat^2):G [Landau diagnostic, ->0]

Mass extraction (massive ansatz D(k) = Z/(khat^2 + m^2)):
    1/D(k) = (1/Z) khat^2 + m^2/Z   -> linear fit, intercept = m^2/Z
    m^2 = intercept/slope, Z = 1/slope, D(0) = 1/intercept = Z/m^2.
    intercept > 0  <=> massive, D(0) finite.   intercept ~ 0 <=> massless.
Gribov-Stingl  D(k) = Z (khat^2+m1^2)/(khat^4+m2^2 khat^2+m3^4).
Jackknife errors over configs.
"""
import sys, time, json
import numpy as np
sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
import w4_blocked as W
from su3 import ID2, SIG
T2 = SIG / 2.0
adj = W.adj


# ---------------------------------------------------------------- gauge fixing
def polar_su2(S):
    Uu, sv, Vh = np.linalg.svd(S)
    return Uu @ Vh


def gauge_fix(U, rng, niter=1500, tol=1e-7):
    L = U.shape[1]
    g = W.random_su2(rng, (L, L, L, L))
    g0 = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
    par = ((g0[0] + g0[1] + g0[2] + g0[3]) % 2).reshape(-1)
    idxp = [np.nonzero(par == c)[0] for c in (0, 1)]
    gflat = g.reshape(-1, 2, 2)

    def Fval(g):
        tot = 0.0
        for mu in range(4):
            tot += np.real(np.trace(g @ U[mu] @ adj(np.roll(g, -1, axis=mu)),
                                    axis1=-2, axis2=-1)).sum()
        return tot

    Fp = None
    for it in range(niter):
        for c in (0, 1):
            B = np.zeros((L, L, L, L, 2, 2), dtype=np.complex128)
            for mu in range(4):
                B += (U[mu] @ adj(np.roll(g, -1, axis=mu))
                      + adj(np.roll(U[mu], 1, axis=mu)) @ adj(np.roll(g, 1, axis=mu)))
            gnew = adj(polar_su2(B))
            gflat[idxp[c]] = gnew.reshape(-1, 2, 2)[idxp[c]]
        F = Fval(g)
        if Fp is not None and abs(F - Fp) < tol * abs(F):
            break
        Fp = F
    out = np.empty_like(U)
    for mu in range(4):
        out[mu] = g @ U[mu] @ adj(np.roll(g, -1, axis=mu))
    return out


def gauge_field(U):
    return np.einsum('...ij,aij->...a', U, T2).imag * 2.0


def div_rms(U):
    A = gauge_field(U)
    div = np.zeros(A.shape[1:], dtype=float)
    for mu in range(4):
        div += A[mu] - np.roll(A[mu], 1, axis=mu)
    return float(np.sqrt(np.mean(div ** 2))), float(np.sqrt(np.mean(A ** 2)))


# ---------------------------------------------------------------- momenta / prop
def kmomenta(L):
    g = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
    khat = [2.0 * np.sin(np.pi * g[mu] / L) for mu in range(4)]
    k2 = sum(khat[mu] ** 2 for mu in range(4))
    return k2, khat


def transverse_prop(U, k2, khat):
    """Returns (D_exact, DL_exact, D_khat, DL_khat).

    The gauge condition imposed by Landau fixing is the lattice divergence
    Sigma_mu [A_mu(x) - A_mu(x-mu)] = 0, whose Fourier form is
        Sigma_mu L_mu At_mu = 0,   L_mu = 1 - exp(-i k_mu).
    The physically transverse projector consistent with THIS condition is
        P^{T,exact}_{mu nu} = delta_{mu nu} - L_mu conj(L_nu)/|L|^2,
    so  D_exact = (1/3)(TrG - <|L.At|^2>/|L|^2),  |L|^2 = khat^2.
    We compute the longitudinal term directly (sum over colours first) which is
    numerically robust; it is ~1e-6 of TrG when gauge fixing has converged.
    D_khat uses the real continuum-style projector delta - khat_mu khat_nu/khat^2
    (the textbook choice) as a cross-check.
    """
    L = U.shape[1]
    V = L ** 4
    A = gauge_field(U)
    At = np.fft.fftn(A, axes=(1, 2, 3, 4))
    Ac = np.conj(At)
    Gmunu = np.einsum('m...c,n...c->mn...', At, Ac) / V
    TrG = sum(Gmunu[m, m] for m in range(4))
    # exact lattice longitudinal vectors L_mu = 1 - exp(-i k_mu)
    n = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
    Lv = [1.0 - np.exp(-2j * np.pi * n[m] / L) for m in range(4)]
    L2 = sum(np.abs(Lv[m]) ** 2 for m in range(4))          # == khat^2
    LA = sum(Lv[m][..., None] * At[m] for m in range(4))    # (L,L,L,L,3)
    lng = (np.abs(LA) ** 2).sum(-1) / V                     # |L.At|^2/V
    safeL = np.where(L2 > 1e-12, L2, 1.0)
    # real-khat projector longitudinal
    kGk = np.zeros(k2.shape, dtype=np.complex128)
    for m in range(4):
        for nn in range(4):
            kGk += khat[m] * khat[nn] * Gmunu[m, nn]
    safeK = np.where(k2 > 1e-12, k2, 1.0)
    D_ex = np.where(k2 > 1e-12, (TrG - lng / safeL) / 3.0, 0.0).real
    DL_ex = np.where(k2 > 1e-12, lng / safeL, 0.0).real
    D_k = np.where(k2 > 1e-12, (TrG - kGk / safeK) / 3.0, 0.0).real
    DL_k = np.where(k2 > 1e-12, kGk / safeK, 0.0).real
    return D_ex, DL_ex, D_k, DL_k


def shell_data(L):
    k2, khat = kmomenta(L)
    mask = k2 > 1e-9
    ushells = np.unique(np.round(k2[mask], 6))
    return k2, khat, ushells


# ---------------------------------------------------------------- free control
def free_config(rng, L, amp=0.02):
    """Sample a FREE, LATTICE-DIVERGENCE-FREE massless field (mimics Landau gauge):
    At(k) lies in the 3-dim complex complement of L_mu=1-exp(-i k_mu), with
    var sigma^2 = amp^2/khat^2 per transverse direction/colour.
    Reconstruct real A, exponentiate. Recover D(k)=amp^2/khat^2 (intercept 0)."""
    k2, khat = kmomenta(L)
    V = L ** 4
    At = np.zeros((4,) + k2.shape + (3,), dtype=np.complex128)
    n = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')

    for k in np.ndindex(*([L] * 4)):
        kk2 = k2[k]
        if kk2 < 1e-9:
            continue
        neg = tuple((L - x) % L for x in k)
        if k > neg:
            continue
        Lv = np.array([1.0 - np.exp(-2j * np.pi * n[m][k] / L) for m in range(4)])
        Vh = np.linalg.svd(Lv.reshape(1, 4))[2]
        e = Vh[1:4]                       # 3 complex orthonormal rows, L-transverse
        sigma = amp / np.sqrt(kk2)
        selfconj = (k == neg)
        for lam in range(3):
            for a in range(3):
                if selfconj:
                    c = sigma * rng.standard_normal()
                else:
                    c = sigma * (rng.standard_normal() + 1j * rng.standard_normal()) / np.sqrt(2.0)
                for m in range(4):
                    At[m][k][a] += c * e[lam][m]
                    if not selfconj:
                        At[m][neg][a] += np.conj(c) * np.conj(e[lam][m])
    A = np.fft.ifftn(At, axes=(1, 2, 3, 4)).real * V
    Af = np.einsum('...a,aij->...ij', A, T2)
    w, Vv = np.linalg.eigh(Af)
    ph = np.exp(1j * w)
    return np.einsum('...ik,...k,...jk->...ij', Vv, ph, np.conj(Vv))


def control_free(L, ncfg, seed=7):
    rng = np.random.default_rng(seed)
    k2, khat, ushells = shell_data(L)
    ns = len(ushells)
    Dcfg = np.full((ncfg, ns), np.nan)
    for c in range(ncfg):
        U = free_config(rng, L)
        D, DL, Dk, DLk = transverse_prop(U, k2, khat)
        for i, sh in enumerate(ushells):
            sel = np.round(k2, 6) == sh
            Dcfg[c, i] = D[sel].mean()
    meanD = np.nanmean(Dcfg, axis=0)
    m2, Z, D0, inter = fit_massive(ushells, meanD)
    print(f"[FREE control L={L}] input massless 1/k2 -> fit: m2={m2:.3e}  Z={Z:.4f}  "
          f"D0={D0:.3e}  intercept={inter:.3e}")
    print("   k2 : D(meas)   [1/k2]")
    for i, sh in enumerate(ushells):
        print(f"   {sh:6.3f} : {meanD[i]:.5f}   [{1/sh:.5f}]")
    return dict(kind='free', L=L, ushells=ushells.tolist(), meanD=meanD.tolist(),
                mass2=m2, Z=Z, D0=D0, intercept=inter)


# ---------------------------------------------------------------- fits
def fit_massive(k2shell, D):
    ok = (D > 1e-12) & np.isfinite(D)
    if ok.sum() < 2:
        return np.nan, np.nan, np.nan, np.nan
    x = k2shell[ok]
    y = 1.0 / D[ok]
    A = np.vstack([x, np.ones_like(x)]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    if abs(a) < 1e-30:
        return np.nan, np.nan, np.nan, np.nan
    m2 = b / a
    Z = 1.0 / a
    D0 = 1.0 / b if b > 1e-30 else np.inf
    return m2, Z, D0, b


def fit_massive_window(k2shell, D, k2max):
    sel = k2shell <= k2max + 1e-9
    return fit_massive(k2shell[sel], D[sel])


def fit_GS(k2shell, D):
    from scipy.optimize import least_squares
    ok = (D > 1e-12) & np.isfinite(D)
    x = k2shell[ok]; y = D[ok]
    if len(x) < 4:
        return np.nan, np.nan
    def res(p):
        Z, m1s, m2s, m3s = p
        return Z * (x + m1s) / (x * x + m2s * x + m3s ** 2) - y
    best = None
    for m1g in (0.05, 0.5, 2.0):
        p0 = [max(y[0] * (x[0] + m1g), 1e-3), m1g, 1.0, 0.5]
        try:
            r = least_squares(res, p0, bounds=(1e-9, np.inf), max_nfev=8000)
            if best is None or r.cost < best.cost:
                best = r
        except Exception:
            pass
    if best is None:
        return np.nan, np.nan
    Z, m1s, m2s, m3s = best.x
    D0 = Z * m1s / (m3s ** 2)
    return D0, np.sqrt(max(m1s, 0.0))


# ---------------------------------------------------------------- chain advance
def _advance(U, beta, rng, n):
    L = U.shape[1]; ns = L ** 4
    g0 = np.meshgrid(*[np.arange(L)] * 4, indexing='ij')
    par = ((g0[0] + g0[1] + g0[2] + g0[3]) % 2).reshape(-1)
    idxp = [np.nonzero(par == c)[0] for c in (0, 1)]
    kappa = beta / 2.0
    Uf = U.reshape(4, ns, 2, 2)
    for _ in range(n):
        for c in (0, 1):
            ids = idxp[c]
            for mu in range(4):
                S = W.staple_sum2(U, mu).reshape(ns, 2, 2)[ids]
                A = Uf[mu, ids]
                v = W.su2_heatbath(kappa, A @ S, 0, 1, rng)
                Uf[mu, ids] = v @ A
    return U


# ---------------------------------------------------------------- driver
def run(beta, L, ncfg, ntherm=500, gap=8, seed=1, k2cut=None):
    rng = np.random.default_rng(seed)
    k2, khat, ushells = shell_data(L)
    ns = len(ushells)
    Dcfg = np.full((ncfg, ns), np.nan)
    DLcfg = np.full((ncfg, ns), np.nan)
    Dkcfg = np.full((ncfg, ns), np.nan)
    divs = np.zeros(ncfg)
    plaqs = np.zeros(ncfg)
    t0 = time.time()
    U = W.sample2(L, beta, rng, ntherm=ntherm, nkeep=1, gapevery=1)[0]
    for c in range(ncfg):
        U = _advance(U, beta, rng, gap)
        Ugf = gauge_fix(U, rng)
        d, rms = div_rms(Ugf)
        divs[c] = d / rms
        plaqs[c] = W.action2(Ugf, beta) / (6.0 * L ** 4)
        Dex, DLex, Dk, DLk = transverse_prop(Ugf, k2, khat)
        for i, sh in enumerate(ushells):
            sel = np.round(k2, 6) == sh
            Dcfg[c, i] = Dex[sel].mean()
            DLcfg[c, i] = DLex[sel].mean()
            Dkcfg[c, i] = Dk[sel].mean()
    meanD = np.nanmean(Dcfg, axis=0)
    meanDL = np.nanmean(DLcfg, axis=0)
    meanDk = np.nanmean(Dkcfg, axis=0)
    if k2cut is None:
        k2cut = ushells[min(len(ushells) - 1, 3)] if len(ushells) > 3 else ushells[-1]

    n = ncfg
    def jk(fn):
        full = np.array(fn(np.arange(n)), dtype=float)
        vals = np.array([np.array(fn(np.delete(np.arange(n), i)), dtype=float)
                         for i in range(n)])
        err = np.sqrt((n - 1) / n * np.sum((vals - vals.mean(axis=0)) ** 2, axis=0))
        return full, err

    def frFull(ids):
        return fit_massive(ushells, np.nanmean(Dcfg[ids], axis=0))
    def frWin(ids):
        return fit_massive_window(ushells, np.nanmean(Dcfg[ids], axis=0), k2cut)
    def frWinK(ids):
        return fit_massive_window(ushells, np.nanmean(Dkcfg[ids], axis=0), k2cut)

    ffull, efull = jk(frFull)
    fwin, ewin = jk(frWin)
    fwink, ewink = jk(frWinK)

    # two-point effective mass from two lowest shells
    def twopt(ids):
        m = np.nanmean(Dcfg[ids], axis=0)
        if m[0] > 1e-12 and m[1] > 1e-12:
            y0, y1 = 1 / m[0], 1 / m[1]
            slope = (y1 - y0) / (ushells[1] - ushells[0])
            inter = y0 - slope * ushells[0]
            return inter / slope if abs(slope) > 1e-30 else np.nan
        return np.nan
    f2, e2 = jk(lambda ids: (twopt(ids),))

    D0g, m1g = fit_GS(ushells, meanD)
    out = dict(kind='mc', beta=beta, L=L, ncfg=ncfg,
               ushells=ushells.tolist(), meanD=meanD.tolist(),
               meanDk=meanDk.tolist(), meanDL=meanDL.tolist(),
               div_rms_ratio=float(np.mean(divs)), plaq=float(np.nanmean(plaqs)),
               k2cut=float(k2cut),
               mass2_full=float(ffull[0]), mass2_full_err=float(efull[0]),
               Z_full=float(ffull[1]), D0_full=float(ffull[2]),
               inter_full=float(ffull[3]),
               mass2_win=float(fwin[0]), mass2_win_err=float(ewin[0]),
               Z_win=float(fwin[1]), D0_win=float(fwin[2]), D0_win_err=float(ewin[2]),
               inter_win=float(fwin[3]),
               mass2_win_khat=float(fwink[0]), mass2_win_khat_err=float(ewink[0]),
               D0_win_khat=float(fwink[2]),
               mass2_2pt=float(f2[0]), mass2_2pt_err=float(e2[0]),
               D0_GS=float(D0g), m1_GS=float(m1g),
               time=time.time() - t0)
    print(f"\n[beta={beta} L={L} ncfg={ncfg}] {out['time']:.1f}s "
          f"div/rms={out['div_rms_ratio']:.2e} plaq={out['plaq']:.4f}")
    print("   k2      D(k)        D_khat      1/D")
    for i, sh in enumerate(ushells):
        print(f"  {sh:6.3f}  {meanD[i]: .5e}  {meanDk[i]: .5e}  {1/meanD[i]: .4f}")
    print(f"   FULL fit : m2={ffull[0]:+.4f}+-{efull[0]:.4f}  m={np.sqrt(max(ffull[0],0)):.4f}"
          f"  Z={ffull[1]:.4f}  D0={ffull[2]:.4f}  intercept={ffull[3]:.4f}")
    print(f"   WIN <= {k2cut:.2f}: m2={fwin[0]:+.4f}+-{ewin[0]:.4f}  m={np.sqrt(max(fwin[0],0)):.4f}"
          f"  Z={fwin[1]:.4f}  D0={fwin[2]:.4f}+-{ewin[2]:.4f}  intercept={fwin[3]:.4f}")
    print(f"   WIN khat: m2={fwink[0]:+.4f}+-{ewink[0]:.4f}  D0={fwink[2]:.4f}")
    print(f"   2pt m2={f2[0]:+.4f}+-{e2[0]:.4f}   GS: D0={D0g:.4f} m1={m1g:.4f}")
    return out


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    results = []
    if which in ('all', 'ctrl'):
        results.append(control_free(4, 20))
        results.append(control_free(6, 12))
    if which in ('all', 'run'):
        for L in (4, 6):
            for beta in (2.2, 2.5, 2.8):
                nc = 120 if L == 4 else int(sys.argv[2]) if len(sys.argv) > 2 else 60
                results.append(run(beta, L, ncfg=nc))
    json.dump(results, open('D:/Dev/math/research/lattice-cpb/gluon_prop_results.json', 'w'), indent=1)
    print("\nWROTE gluon_prop_results.json")
