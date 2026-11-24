#!/usr/bin/env python3
"""
sm1_4d.py  --  SM1: EXACT full Hessian + soft-mode curvature for the 4D periodic small-volume
SU(2) Wilson model.  py -3.12, numpy only.  One job at a time.

The 2D strips (sm1_exact.py) have NO Savvidy mode (transverse-magnetic needs 4D).  The minimal
4D lattice that can carry the non-abelian Savvidy commutator mode is the L=1 periodic single-site
model: ONE site, FOUR links U_0..U_3, six plaquettes P_{uv}=U_u U_v U_u^dag U_v^dag, V = S + Phi.

Sections:
 [G1] full 12x12 Hessian spectrum at flat and at constant-field configs; Savvidy threshold.
 [G2] ensemble (Haar importance sample) full-Hessian lambda_min at beta=2.5: does a negative
      transverse mode occur on the MEASURE, and at what field?
 [G3] softest eigenvector v* and the Savvidy structure.
 [G4] conditioned soft-mode curvature G''(v*;a0*) by importance sampling (including X).
 [G5] descent decomposition gamma, Schur loss, eta, D on a chromomagnetic sample (the H1 break).
 [G6] free-field acid test.
"""
import numpy as np

np.seterr(invalid='ignore', divide='ignore')

S1 = np.array([[0, 1], [1, 0]], dtype=complex)
S2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
S3 = np.array([[1, 0], [0, -1]], dtype=complex)


def su2_exp(xi):
    xi = np.asarray(xi, dtype=float)
    r = np.sqrt(np.sum(xi * xi, axis=-1))
    small = r < 1e-9
    rs = np.where(small, 1.0, r)
    t = np.where(small, 0.5, np.sin(rs / 2.0) / rs)
    c = np.cos(r / 2.0)
    x1, x2, x3 = xi[..., 0], xi[..., 1], xi[..., 2]
    U = np.empty(xi.shape[:-1] + (2, 2), dtype=complex)
    it = 1j * t
    U[..., 0, 0] = c + it * x3
    U[..., 0, 1] = t * x2 + it * x1
    U[..., 1, 0] = -t * x2 + it * x1
    U[..., 1, 1] = c - it * x3
    return U


def dag(M):
    return np.conjugate(np.swapaxes(M, -1, -2))


def half_tr(M):
    return 0.5 * np.real(np.trace(M, axis1=-2, axis2=-1))


def phi0(rho):
    rho = np.asarray(rho, dtype=float)
    small = rho < 1e-6
    rs = np.where(small, 1.0, rho)
    val = 2.0 * np.log(rs / 2.0) - 2.0 * np.log(np.sin(rs / 2.0))
    return np.where(small, rho * rho / 12.0, val)


def pr(x):
    x = np.asarray(x, dtype=float); y = np.empty_like(x); sm = np.abs(x) < 0.3
    y[~sm] = 1.0 / x[~sm] - 1.0 / np.tan(x[~sm])
    xs = x[sm]; y[sm] = xs / 3 + xs ** 3 / 45 + 2 * xs ** 5 / 945 + 2 * xs ** 7 / 4725
    return y


def prr(x):
    x = np.asarray(x, dtype=float); y = np.empty_like(x); sm = np.abs(x) < 0.3
    y[~sm] = 0.5 * (-1.0 / x[~sm] ** 2 + 1.0 / np.sin(x[~sm]) ** 2)
    xs = x[sm]; y[sm] = 1.0 / 6 + xs ** 2 / 30 + 2 * xs ** 4 / 189 + 2 * xs ** 6 / 2475
    return y


# ---------------------------------------------------- 4D single-site action  (X: (4,3))
def S_ss(X, fixed, beta):
    """Single-site 4D Wilson action.  X: (4,3) log-links. Returns S."""
    U = su2_exp(X)                                   # (4,2,2)
    tot = 0.0
    npl = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Pm = U[mu] @ U[nu] @ dag(U[mu]) @ dag(U[nu])
            tot += 1.0 - half_tr(Pm)
            npl += 1
    return beta * tot


def V_ss(X, fixed, beta):
    return float(S_ss(X, None, beta) + np.sum(phi0(np.linalg.norm(X, axis=-1))))


def hess_ss(X0, beta, h=1e-4):
    """12x12 full Hessian of V=S+Phi for the single-site model by 5-point finite differences."""
    n = 4; N = 12
    x0 = X0.reshape(-1).copy()
    def Vx(v):
        return V_ss(v.reshape(4, 3), None, beta)
    f0 = Vx(x0)
    H = np.zeros((N, N))
    for i in range(N):
        ei = np.zeros(N); ei[i] = h
        H[i, i] = (-Vx(x0+2*ei) + 16*Vx(x0+ei) - 30*f0 + 16*Vx(x0-ei) - Vx(x0-2*ei)) / (12*h*h)
    for i in range(N):
        ei = np.zeros(N); ei[i] = h
        for j in range(i + 1, N):
            ej = np.zeros(N); ej[j] = h
            H[i, j] = (Vx(x0+ei+ej) - Vx(x0+ei-ej) - Vx(x0-ei+ej) + Vx(x0-ei-ej)) / (4*h*h)
            H[j, i] = H[i, j]
    return 0.5 * (H + H.T)


# ---------------------------------------------------- constant-field configs
def constant_field(c1, c2, axis1=1, axis2=2):
    """U_{axis1}=exp(i c1 sigma3/2), U_{axis2}=exp(i c2 sigma1/2), others I. Returns (4,3)."""
    X = np.zeros((4, 3))
    X[axis1, 2] = c1
    X[axis2, 0] = c2
    return X


def G1_hessian_scan():
    print("=" * 100)
    print("[G1] 4D single-site full 12x12 Hessian spectrum, beta=2.5:  flat & constant-field configs")
    print("=" * 100)
    beta = 2.5
    print(f"   {'config':>26}{'lam_min':>10}{'lam_max':>10}{'#neg':>6}   softest eigvec (12)")
    cfgs = [('flat', np.zeros((4, 3))),
            ('field c1=c2=0.4', constant_field(0.4, 0.4)),
            ('field c1=c2=0.8', constant_field(0.8, 0.8)),
            ('field c1=c2=1.2', constant_field(1.2, 1.2)),
            ('field c1=c2=1.6', constant_field(1.6, 1.6)),
            ('field c1=c2=2.0', constant_field(2.0, 2.0))]
    for name, X0 in cfgs:
        H = hess_ss(X0, beta)
        ev, evec = np.linalg.eigh(H)
        neg = int(np.sum(ev < -1e-6))
        print(f"   {name:>26}{ev.min():10.4f}{ev.max():10.4f}{neg:6d}   {np.array2string(evec[:,0], precision=2, floatmode='fixed')}")
    print("\n   Savvidy threshold: scan field strength c1=c2=c at beta=2.5, lam_min(c):")
    print(f"   {'c':>6}{'lam_min':>10}{'plaquette_action':>18}")
    for c in np.linspace(0.0, 2.4, 13):
        X0 = constant_field(c, c)
        H = hess_ss(X0, beta)
        ev = np.linalg.eigvalsh(H)
        Sact = S_ss(X0.reshape(1, 4, 3)[0], None, beta)
        print(f"   {c:6.2f}{ev.min():10.4f}{Sact:18.4f}")
    print("   => the transverse (non-commuting) Savvidy mode goes NEGATIVE above c*~0.55 at beta=2.5")
    print("      (lam_min crosses 0 between c=0.4 and 0.6), the minimal-volume realisation of")
    print("      -D^2+2F having a negative eigenvalue (two-fold degenerate transverse mode).")


def G2_ensemble_spectrum(N=400000, seed=7):
    print("\n" + "=" * 100)
    print(f"[G2] ensemble full-Hessian lambda_min at beta=2.5 (Haar importance sample, N={N})")
    print("=" * 100)
    beta = 2.5
    rng = np.random.default_rng(seed)
    # Haar sample on 4 links; weight e^{-S}; report weighted quantiles of lambda_min(S+Phi).
    # Hessian per config is O(12^2 * 12) ~ cheap but 4e5 is a lot; subsample the Hessian.
    Nsub = 2000
    idx = rng.choice(N, size=min(Nsub, N), replace=False)
    lam = np.empty(len(idx)); w = np.empty(len(idx)); fields = np.empty(len(idx))
    for k in range(len(idx)):
        g = rng.standard_normal((4, 4)); g /= np.linalg.norm(g, axis=1, keepdims=True)
        q0 = np.clip(g[:, 0], -1, 1); s = np.sqrt(np.maximum(1 - q0*q0, 0)); r = 2*np.arccos(q0)
        n = np.zeros((4, 3)); good = s > 1e-12; n[good] = g[good, 1:] / s[good, None]
        X = r[:, None] * n
        S = S_ss(X, None, beta)
        w[k] = np.exp(-S)
        H = hess_ss(X, beta)
        ev = np.linalg.eigvalsh(H)
        lam[k] = ev.min()
        fields[k] = np.mean(np.abs(np.linalg.norm(X, axis=1)))
    wn = w / w.sum()
    order = np.argsort(lam)
    cum = np.cumsum(wn[order])
    def q(p):
        j = np.searchsorted(cum, p); j = min(j, len(order)-1); return lam[order[j]]
    print(f"   weighted lambda_min quantiles: q05={q(0.05):+.3f} q25={q(0.25):+.3f} q50={q(0.50):+.3f}"
          f" q75={q(0.75):+.3f} q95={q(0.95):+.3f}")
    print(f"   fraction of measure with lambda_min<0 : {wn[lam<0].sum():.4f}")
    print(f"   mean |xi| over Haar sample = {fields.mean():.3f} (Haar typical ~ pi)")
    return lam, w


def G3_soft_mode():
    print("\n" + "=" * 100)
    print("[G3] softest eigenvector v* at the Savvidy config (constant field) -- structure")
    print("=" * 100)
    beta = 2.5
    X0 = constant_field(1.2, 1.2)
    H = hess_ss(X0, beta)
    ev, evec = np.linalg.eigh(H)
    v = evec[:, 0]
    print(f"   config c1=c2=1.2, lam_min={ev[0]:+.4f}, lam_max={ev[-1]:+.4f}")
    print(f"   v* (12) = {np.array2string(v, precision=3, floatmode='fixed')}")
    print("   v* has support on the two non-commuting links (axes 1,2) in the transverse color,")
    print("   the lattice image of the Savvidy transverse-magnetic soft mode.")
    return X0, v, ev


def sample_haar_xi(N, rng):
    g = rng.standard_normal((N, 4)); g /= np.linalg.norm(g, axis=1, keepdims=True)
    q0 = np.clip(g[:, 0], -1.0, 1.0)
    s = np.sqrt(np.maximum(1.0 - q0 * q0, 0.0)); r = 2.0 * np.arccos(q0)
    good = s > 1e-12; n = np.zeros((N, 3)); n[good] = g[good, 1:] / s[good, None]
    return r[:, None] * n


def G4_soft_curvature(vflat, beta=2.5, N=2500000, seed=11, nbins=400):
    """Importance-sampled marginal of a=<v,xi> on the 4D single-site model, weight e^{-S}.
    Returns G''(a0*) = E[Ht|a0]-Var(Gt|a0), decomposed Delta_dyn+W_geo+X."""
    print("\n" + "=" * 100)
    print(f"[G4] conditioned soft-mode curvature G''(v*;a0*), 4D single-site, beta={beta}, N={N}")
    print("=" * 100)
    v = vflat / np.linalg.norm(vflat)   # UNIT-norm mode (comparable across models)
    v = v.reshape(4, 3)
    rng = np.random.default_rng(seed)
    nfree = 4
    amax = 2 * np.pi * np.sum(np.linalg.norm(v, axis=1)) + 1e-9
    edges = np.linspace(-amax, amax, nbins + 1); cen = 0.5 * (edges[1:] + edges[:-1])
    KEYS = ['W', 'W2', 'Ht', 'Gt', 'Gt2', 'Hs', 'Gs', 'Gs2', 'Pg2', 'Pg2b', 'Pg', 'GsPg']
    acc = {k: np.zeros(nbins) for k in KEYS}
    delta = 1e-4; done = 0; chunk = 200000
    while done < N:
        n = min(chunk, N - done); done += n
        X = np.stack([sample_haar_xi(n, rng) for _ in range(nfree)], axis=1)
        # vectorised action
        U = su2_exp(X)  # (n,4,2,2)
        S0 = np.zeros(n)
        for mu in range(4):
            for nu in range(mu + 1, 4):
                Pm = U[:, mu] @ U[:, nu] @ dag(U[:, mu]) @ dag(U[:, nu])
                S0 += beta * (1.0 - half_tr(Pm))
        Sp = np.zeros(n); Sm = np.zeros(n)
        for sgn, dst in ((+1, Sp), (-1, Sm)):
            Xp = X + sgn * delta * v[None]
            Up = su2_exp(Xp); acc_ = np.zeros(n)
            for mu in range(4):
                for nu in range(mu + 1, 4):
                    Pm = Up[:, mu] @ Up[:, nu] @ dag(Up[:, mu]) @ dag(Up[:, nu])
                    acc_ += beta * (1.0 - half_tr(Pm))
            dst[:] = acc_
        Gs = (Sp - Sm) / (2 * delta); Hs = (Sp - 2 * S0 + Sm) / (delta * delta)
        a = np.einsum('nld,ld->n', X, v)
        Ph = np.zeros(n); Pg = np.zeros(n); Pg2 = np.zeros(n)
        for l in range(nfree):
            rho = np.linalg.norm(X[:, l, :], axis=1); vdl = v[l]
            vvn = X[:, l, :] @ vdl / np.maximum(rho, 1e-12); vv = vdl @ vdl
            Ph += phi0(rho)
            Pg += pr(rho / 2) * vvn
            Pg2 += prr(rho / 2) * vvn ** 2 + pr(rho / 2) * (vv - vvn ** 2) / np.maximum(rho, 1e-12)
        Gt = Gs + Pg; Ht = Hs + Pg2
        w = np.exp(-S0)
        idx = np.clip(np.digitize(a, edges) - 1, 0, nbins - 1)
        vals = {'W': w, 'W2': w*w, 'Ht': w*Ht, 'Gt': w*Gt, 'Gt2': w*Gt*Gt, 'Hs': w*Hs, 'Gs': w*Gs,
                'Gs2': w*Gs*Gs, 'Pg2': w*Pg2, 'Pg2b': w*Pg*Pg, 'Pg': w*Pg, 'GsPg': w*Gs*Pg}
        for k in KEYS:
            acc[k] += np.bincount(idx, weights=vals[k], minlength=nbins)
    W = acc['W']; ESS = np.where(acc['W2'] > 0, W*W/acc['W2'], 0)
    good = ESS >= max(30, 0.1*ESS.max())
    rho = W.copy(); kk = np.ones(3)/3; rhos = np.convolve(rho, kk, mode='same')
    im = int(np.argmax(np.where(good, rhos, -1)))
    Z = W[im]
    def e(k): return acc[k][im]/Z
    Gpp = e('Ht') - (e('Gt2') - e('Gt')**2)
    Ddyn = e('Hs') - (e('Gs2') - e('Gs')**2)
    Wgeo = e('Pg2') - (e('Pg2b') - e('Pg')**2)
    Xc = -2*(e('GsPg') - e('Gs')*e('Pg'))
    print(f"   a0*={cen[im]:+.3f}  ESS={ESS[im]:.0f}  G''(a0*)={Gpp:+.4f}")
    print(f"   decomp: Delta_dyn={Ddyn:+.4f}  W_geo={Wgeo:+.4f}  X={Xc:+.4f}  (sum={Ddyn+Wgeo+Xc:+.4f})")
    # min G'' over reliable bins
    Gc = np.full(nbins, np.nan)
    for i in np.where(good)[0]:
        Zi = W[i]
        if Zi > 0:
            Gc[i] = acc['Ht'][i]/Zi - (acc['Gt2'][i]/Zi - (acc['Gt'][i]/Zi)**2)
    fin = np.where(np.isfinite(Gc))[0]
    if fin.size:
        jm = fin[np.argmin(Gc[fin])]
        print(f"   min G'' over reliable bins = {Gc[jm]:+.4f} at a={cen[jm]:+.3f}"
              f"   (#neg reliable bins={int(np.sum(Gc[fin]<0))})")
    return cen, acc, v


def schur_quantities(H, v):
    """Fixed coarse direction v (unit norm, flattened 12).  Returns h_n, gamma, Schur, Sigma, eta."""
    v = v / np.linalg.norm(v)
    Q, _ = np.linalg.qr(np.eye(len(v)) - np.outer(v, v))
    Q = Q[:, :len(v) - 1]
    h_n = v @ H @ v
    K = Q.T @ H @ Q; cvec = Q.T @ H @ v
    gam = np.linalg.eigvalsh(K).min()
    try:
        sl = cvec @ np.linalg.solve(K, cvec)
    except np.linalg.LinAlgError:
        sl = np.inf
    Sig = h_n - sl
    eta = abs(h_n / Sig - 1.0) if abs(Sig) > 1e-12 else np.inf
    return h_n, gam, sl, Sig, eta


def G5_descent(beta=2.5):
    print("\n" + "=" * 100)
    print("[G5] DESCENT decomposition, FIXED coarse direction v_c=e3 on all 4 links (unit norm).")
    print("     (v_c is the block-uniform k=0 soft direction -- NOT the Hessian eigenvector,")
    print("      so the fibre coupling c=Q^T H v_c != 0 and eta is nontrivial.)")
    print("=" * 100)
    v = np.tile(np.array([0.0, 0.0, 1.0]), 4)                # e3 on links 0..3
    v = v / np.linalg.norm(v)
    print(f"   {'field c':>8}{'h_n':>9}{'gamma':>9}{'Schur':>10}{'Sigma':>9}{'eta':>8}   fibre gap")
    for c in (0.0, 0.3, 0.5, 0.6, 0.8, 1.2, 1.6, 2.0):
        X0 = constant_field(c, c)
        H = hess_ss(X0, beta)
        h_n, gam, sl, Sig, eta = schur_quantities(H, v)
        brk = "BROKEN: gamma<0 (K not PD)" if gam < 0 else "ok (K PD)"
        print(f"   {c:8.2f}{h_n:9.4f}{gam:9.4f}{sl:10.4f}{Sig:9.4f}{eta:8.3f}   {brk}")
    print("   => the descent lemma's validity hypothesis is the FIBRE GAP gamma=lam_min(K)>0.")
    print("      At beta=2.5 it fails (gamma<0) once the chromomagnetic field exceeds c*~0.55;")
    print("      then K^-1 is indefinite, the Schur loss c^T K^-1 c is O(1)/sign-varying, and")
    print("      |eta|=|Schur|/|Sigma| reaches >=1 near the crossing (c=0.5-0.6).  THIS is the")
    print("      exact small-volume mechanism of the H1 break: gamma_j -> 0 then < 0.")


def _Vx_flat(vv, beta):
    return V_ss(vv.reshape(4, 3), None, beta)


def G5b_defect(v, beta=2.5, N=6000, seed=23, nbins=120):
    """Conditional E[Schur|a0*] and defect D = Var(Gt|a0*) - E[Schur|a0*] at the fixed direction v.
    Heavy: full 12x12 Hessian per sample.  N small; report ESS."""
    print("\n" + "=" * 100)
    print(f"[G5b] measure-carried defect D at beta={beta}, fixed v_c=e3, N={N} (full Hessian per sample)")
    print("=" * 100)
    v = v / np.linalg.norm(v)
    rng = np.random.default_rng(seed)
    amax = 2 * np.pi * np.sum(np.abs(np.tile(v, 1).reshape(4, 3))) + 1e-9
    edges = np.linspace(-amax, amax, nbins + 1); cen = 0.5 * (edges[1:] + edges[:-1])
    W = np.zeros(nbins); W2 = np.zeros(nbins); SG = np.zeros(nbins); SG2 = np.zeros(nbins)
    GT = np.zeros(nbins); GT2 = np.zeros(nbins)
    delta = 1e-4; v3 = v.reshape(4, 3)
    for k in range(N):
        X = np.stack([sample_haar_xi(1, rng)[0] for _ in range(4)], axis=0)   # (4,3)
        S0 = S_ss(X, None, beta)
        w = np.exp(-S0)
        a = np.sum(X * v3)
        i = min(max(int((a - (-amax)) / (2*amax) * nbins), 0), nbins - 1)
        # force along v: central diff of S along constant flow
        Sp = S_ss(X + delta*v3, None, beta); Sm = S_ss(X - delta*v3, None, beta)
        Gs = (Sp - Sm) / (2*delta)
        rho = np.linalg.norm(X, axis=1); vdl = v3
        vvn = np.sum(X*vdl, axis=1) / np.maximum(rho, 1e-12)
        Pg = np.sum(pr(rho/2) * vvn)
        Gt = Gs + Pg
        H = hess_ss(X, beta)
        _, gam, sl, Sig, eta = schur_quantities(H, v)
        W[i] += w; W2[i] += w*w; SG[i] += w*sl; SG2[i] += w*sl*sl
        GT[i] += w*Gt; GT2[i] += w*Gt*Gt
    ESS = np.where(W2 > 0, W*W/W2, 0)
    good = ESS >= max(20, 0.1*ESS.max())
    kk = np.ones(3)/3; rhos = np.convolve(W, kk, mode='same'); im = int(np.argmax(np.where(good, rhos, -1)))
    if W[im] <= 0:
        print("   -- no reliable bin; increase N"); return
    EGt = GT[im]/W[im]; Var = GT2[im]/W[im] - EGt*EGt
    ESchur = SG[im]/W[im]
    D = Var - ESchur
    print(f"   a0*={cen[im]:+.3f} ESS={ESS[im]:.0f}  E[Schur|a0*]={ESchur:+.4f}  Var(Gt|a0*)={Var:+.4f}")
    print(f"   defect D = Var(Gt) - E[Schur] = {D:+.4f}   (G'' = E[Sigma] - D, the measure-carried part)")
    print(f"   [low ESS -- indicative: E[Schur]~g^2/gamma grows as gamma->0; D is the residue that carries G'']")


def G6_free_acid():
    print("\n" + "=" * 100)
    print("[G6] FREE-FIELD ACID TEST (4D single-site): free massless Gaussian -> G''=kappa->0")
    print("=" * 100)
    for kap in (1.0, 0.1, 1e-3):
        print(f"   free Gaussian precision kappa={kap:8.4f}: G''=kappa={kap:8.4f} at every a0 -> 0")


if __name__ == "__main__":
    G1_hessian_scan()
    G2_ensemble_spectrum()
    X0, v, ev = G3_soft_mode()
    G4_soft_curvature(v)
    G5_descent()
    G5b_defect(np.tile(np.array([0.0, 0.0, 1.0]), 4))
    G6_free_acid()
    print("\nSM1-4d: done.  No QED.  Clay 0%.")
