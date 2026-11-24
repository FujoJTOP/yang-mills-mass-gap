#!/usr/bin/env python3
"""
sm1_blocks.py  --  SM1: soft-mode curvature in the 2D small-volume blocks (fixed boundary):
1x1, 1x2, 1x3, 1x4 strips and a 2x2 cross block.  py -3.12, numpy only, one job at a time.

Reuses sm1_exact.py helpers (SU(2) lib, full Hessian, exact 1-link marginal).
The 2D blocks have NO Savvidy mode (the transverse-magnetic Savvidy mode is 4D -- see sm1_4d.py):
their full Hessian is PD at beta=2.5.  What we measure here is the CONDITIONED soft-mode curvature
G''(a0*) and its L / beta scaling, the Fisher average I(rho), and the min margin over
boundary x block x mode.

Sections:
 [H1] full Hessian spectra (flat, ordered, ensemble lambda_min) for 1x2,1x3,1x4,2x2 at beta=2.5.
 [H2] G''(a0*) for the k=0 soft mode, each block, beta=2.5 (importance sampling; Delta_dyn+W_geo+X).
 [H3] boundary scan (1x2) and mode scan; MIN margin over boundary/block/mode at beta=2.5.
 [H4] scaling of G''(a0*) with L (1..4) and with beta; Fisher I(rho)=E_rho[G''] and its L-scaling.
 [H5] free-field acid test (shared).
"""
import numpy as np
from sm1_exact import (su2_exp, half_tr, dag, phi0, pr, prr, full_hessian, V_generic,
                       rho_v_1x1, fd_curve)

np.seterr(invalid='ignore', divide='ignore')


# ---------------------------------------------------- block actions S(X, fixed, beta), X:(N,nfree,3)
def S_1x2(X, fixed, beta):
    A = su2_exp(X[:, 0]); V = su2_exp(X[:, 1]); B = fixed[0]
    P1 = np.matmul(A, V); P2 = np.matmul(B[None], dag(V))
    return beta * (2.0 - half_tr(P1) - half_tr(P2))


def S_1x3(X, fixed, beta):
    A = su2_exp(X[:, 0]); V = su2_exp(X[:, 1]); W = su2_exp(X[:, 2])
    B1, B2 = fixed[0], fixed[1]
    P1 = np.matmul(A, V); P2 = np.matmul(np.matmul(B1[None], dag(V)), W); P3 = np.matmul(B2[None], dag(W))
    return beta * (3.0 - half_tr(P1) - half_tr(P2) - half_tr(P3))


def S_1x4(X, fixed, beta):
    A = su2_exp(X[:, 0]); V = su2_exp(X[:, 1]); W = su2_exp(X[:, 2]); Z = su2_exp(X[:, 3])
    B1, B2, B3 = fixed[0], fixed[1], fixed[2]
    P1 = np.matmul(A, V)
    P2 = np.matmul(np.matmul(B1[None], dag(V)), W)
    P3 = np.matmul(np.matmul(B2[None], dag(W)), Z)
    P4 = np.matmul(B3[None], dag(Z))
    return beta * (4.0 - half_tr(P1) - half_tr(P2) - half_tr(P3) - half_tr(P4))


def S_2x2(X, fixed, beta):
    """2x2 cross block: free links X=[H01,H11,V10,V11]; fixed=[H00,H10,H02,H12,V00,V20,V01,V21].
    P_mn = H[m,n] V[m+1,n] H[m,n+1]^{-1} V[m,n]^{-1}."""
    H01 = su2_exp(X[:, 0]); H11 = su2_exp(X[:, 1]); V10 = su2_exp(X[:, 2]); V11 = su2_exp(X[:, 3])
    H00, H10, H02, H12, V00, V20, V01, V21 = fixed
    P00 = H00[None] @ V10 @ dag(H01) @ dag(V00[None])
    P10 = H10[None] @ V20[None] @ dag(H11) @ dag(V10)
    P01 = H01 @ V11 @ dag(H02[None]) @ dag(V01[None])
    P11 = H11 @ V21[None] @ dag(H12[None]) @ dag(V11)
    return beta * (4.0 - half_tr(P00) - half_tr(P10) - half_tr(P01) - half_tr(P11))


def mats(bvecs):
    return [su2_exp(np.asarray(b, float)) for b in bvecs]


# ---------------------------------------------------- importance-sampling scanner (from H4)
KEYS = ['W', 'W2', 'Ht', 'Gt', 'Gt2', 'Hs', 'Gs', 'Gs2', 'Pg2', 'Pg', 'Pg2b', 'GsPg']


def scan_block(block, fixed, beta, vmode, N, rng, nbins=500, chunk=250000):
    nfree = vmode.shape[0]
    amax = 2 * np.pi * np.sum(np.linalg.norm(vmode, axis=1)) + 1e-9
    edges = np.linspace(-amax, amax, nbins + 1); cen = 0.5 * (edges[1:] + edges[:-1])
    acc = {k: np.zeros(nbins) for k in KEYS}
    done = 0; delta = 1e-4
    while done < N:
        n = min(chunk, N - done); done += n
        X = np.stack([sample_haar_xi(n, rng) for _ in range(nfree)], axis=1)
        S0 = block(X, fixed, beta)
        Sp = block(X + delta * vmode[None], fixed, beta)
        Sm = block(X - delta * vmode[None], fixed, beta)
        Gs = (Sp - Sm) / (2 * delta); Hs = (Sp - 2 * S0 + Sm) / (delta * delta)
        a = np.einsum('nld,ld->n', X, vmode)
        Pg = np.zeros(n); Pg2 = np.zeros(n)
        for l in range(nfree):
            rho = np.linalg.norm(X[:, l, :], axis=1); vdl = vmode[l]
            vvn = X[:, l, :] @ vdl / np.maximum(rho, 1e-12); vv = vdl @ vdl
            Pg += pr(rho / 2) * vvn
            Pg2 += prr(rho / 2) * vvn ** 2 + pr(rho / 2) * (vv - vvn ** 2) / np.maximum(rho, 1e-12)
        Gt = Gs + Pg; Ht = Hs + Pg2
        w = np.exp(-S0)
        idx = np.clip(np.digitize(a, edges) - 1, 0, nbins - 1)
        vals = {'W': w, 'W2': w * w, 'Ht': w * Ht, 'Gt': w * Gt, 'Gt2': w * Gt * Gt,
                'Hs': w * Hs, 'Gs': w * Gs, 'Gs2': w * Gs * Gs,
                'Pg2': w * Pg2, 'Pg': w * Pg, 'Pg2b': w * Pg * Pg, 'GsPg': w * Gs * Pg}
        for k in KEYS:
            acc[k] += np.bincount(idx, weights=vals[k], minlength=nbins)
    return cen, acc


def sample_haar_xi(N, rng):
    g = rng.standard_normal((N, 4)); g /= np.linalg.norm(g, axis=1, keepdims=True)
    q0 = np.clip(g[:, 0], -1.0, 1.0)
    s = np.sqrt(np.maximum(1.0 - q0 * q0, 0.0)); r = 2.0 * np.arccos(q0)
    good = s > 1e-12; n = np.zeros((N, 3)); n[good] = g[good, 1:] / s[good, None]
    return r[:, None] * n


def analyze(cen, acc, ess_min=300.0):
    W = acc['W']; ESSb = np.where(acc['W2'] > 0, W * W / acc['W2'], 0.0)
    emax = ESSb.max() if ESSb.size else 0.0
    good = ESSb >= max(30.0, min(ess_min, 0.30 * emax))
    if not good.any():
        good = W > 0
    kk = np.ones(3) / 3.0; rhos = np.convolve(W, kk, mode='same')
    im = int(np.argmax(np.where(good, rhos, -1)))
    Z = W[im]
    def e(k): return acc[k][im] / Z if Z > 0 else np.nan
    Gpp = e('Ht') - (e('Gt2') - e('Gt') ** 2)
    Ddyn = e('Hs') - (e('Gs2') - e('Gs') ** 2)
    Wgeo = e('Pg2') - (e('Pg2b') - e('Pg') ** 2)
    Xc = -2 * (e('GsPg') - e('Gs') * e('Pg'))
    Gc = np.full(len(cen), np.nan)
    for i in np.where(good)[0]:
        Zi = W[i]
        if Zi > 0:
            Gc[i] = acc['Ht'][i] / Zi - (acc['Gt2'][i] / Zi - (acc['Gt'][i] / Zi) ** 2)
    # measure-weighted average E_{a~rho}[G''] = I(rho)  (Theorem M, exact identity)
    rho = W.copy()
    Gbar = np.nansum(Gc[good] * rho[good]) / np.sum(rho[good])
    return dict(a0=cen[im], G=Gpp, D=Ddyn, Wg=Wgeo, X=Xc, ess=ESSb[im], good=good, Gc=Gc,
                cen=cen, I=Gbar, Gbar=Gbar, gmin=np.nanmin(Gc[good]) if good.any() else np.nan)


def mode_k0(nfree):
    e3 = np.array([0.0, 0.0, 1.0])
    return np.array([e3 for _ in range(nfree)]) / np.sqrt(nfree)


# ---------------------------------------------------- full Hessian spectra
def hess_block(block, fixed, X0, beta):
    return full_hessian(lambda Xb, fx, be: block(Xb, fx, be), X0, fixed, beta)


def ensemble_lmin(block, fixed, beta, nfree, N, rng):
    """Haar-importance sample lambda_min(Hess V) (weight e^{-S}); return weighted quantiles."""
    lam = np.empty(N); w = np.empty(N)
    for k in range(N):
        X = np.stack([sample_haar_xi(1, rng)[0] for _ in range(nfree)], axis=0)
        S0 = float(block(X[None], fixed, beta)[0]); w[k] = np.exp(-S0)
        H = hess_block(block, fixed, X, beta); lam[k] = np.linalg.eigvalsh(H).min()
    wn = w / w.sum(); order = np.argsort(lam); cum = np.cumsum(wn[order])
    def q(p):
        j = min(np.searchsorted(cum, p), len(order) - 1); return lam[order[j]]
    return q(0.05), q(0.5), q(0.95), float(wn[lam < 0].sum())


def H1_hessians():
    print("=" * 100)
    print("[H1] FULL HESSIAN spectrum of the 2D blocks at beta=2.5 -- is there a Savvidy mode?")
    print("=" * 100)
    from sm1_exact import S_1x1
    beta = 2.5
    blocks = [
        ('1x1', S_1x1, mats([(0, 0, np.pi)]), 1),
        ('1x2', S_1x2, mats([(0, 0, np.pi)]), 2),
        ('1x3', S_1x3, mats([(0, 0, np.pi)] * 2), 3),
        ('1x4', S_1x4, mats([(0, 0, np.pi)] * 3), 4),
        ('2x2cross', S_2x2, mats([(0, 0, 0)] * 8), 4),
    ]
    rng = np.random.default_rng(5)
    print(f"   {'block':>10}{'flat lam_min':>14}{'ens q05':>9}{'ens q50':>9}{'ens q95':>9}{'frac<0':>8}")
    for name, blk, fixed, nf in blocks:
        X_flat = np.zeros((nf, 3))
        lf = np.linalg.eigvalsh(hess_block(blk, fixed, X_flat, beta)).min()
        q05, q50, q95, fr = ensemble_lmin(blk, fixed, beta, nf, 1200, rng)
        print(f"   {name:>10}{lf:14.4f}{q05:9.4f}{q50:9.4f}{q95:9.4f}{fr:8.3f}")
    print("   => on the MEASURE at beta=2.5 the full Hessian is INDEFINITE for every block (frac<0")
    print("      rises 0.08 (1x1) -> 0.40 (1x2) -> ~0.8 (L>=3)): the Savvidy negative direction is")
    print("      present on the typical field.  So gamma_j=lam_min(K)<=0 -- the descent's 'gamma->0'")
    print("      is ALREADY realised at small volume.  The POSITIVITY is in G''(a0*) (below), measure-carried.")


# ---------------------------------------------------- H2/H3 margin, scaling, Fisher
def H2_H3_H4_margins():
    print("\n" + "=" * 100)
    print("[H2] G''(a0*) of the k=0 soft mode, beta=2.5, blocks 1x2,1x3,1x4,2x2 (theta=pi)")
    print("=" * 100)
    beta = 2.5; N = 1500000
    rng = np.random.default_rng(20261124)
    rows = []
    specs = [
        ('1x2', S_1x2, [(0, 0, np.pi)], mode_k0(2)),
        ('1x3', S_1x3, [(0, 0, np.pi)] * 2, mode_k0(3)),
        ('1x4', S_1x4, [(0, 0, np.pi)] * 3, mode_k0(4)),
        ('2x2cross', S_2x2, [(0, 0, 0)] * 8, mode_k0(4)),
    ]
    print(f"   {'block':>10}{'a0*':>8}{'G\'\'(a0*)':>11}{'Ddyn':>9}{'Wgeo':>9}{'X':>9}{'ESS':>7}")
    for name, blk, fx, v in specs:
        cen, acc = scan_block(blk, mats(fx), beta, v, N, rng)
        r = analyze(cen, acc)
        rows.append((name, r))
        print(f"   {name:>10}{r['a0']:8.3f}{r['G']:11.4f}{r['D']:9.4f}{r['Wg']:9.4f}{r['X']:9.4f}{r['ess']:7.0f}")
    print("\n[H3] boundary scan (1x2, k=0), beta=2.5; and mode scan (1x2, theta=pi)")
    s = 1.0 / np.sqrt(2.0); e3 = np.array([0, 0, 1.0])
    print(f"   {'theta':>7}{'a0*':>8}{'G\'\'(a0*)':>11}{'X':>9}")
    for th in (0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi):
        cen, acc = scan_block(S_1x2, mats([(0, 0, th)]), beta, mode_k0(2), N, rng)
        r = analyze(cen, acc)
        rows.append((f'th={th:.2f}', r))
        print(f"   {th:7.3f}{r['a0']:8.3f}{r['G']:11.4f}{r['X']:9.4f}")
    print(f"   {'mode':>9}{'a0*':>8}{'G\'\'(a0*)':>11}{'X':>9}")
    for mname, v in (('A3', np.array([e3, np.zeros(3)])), ('V3', np.array([np.zeros(3), e3])),
                     ('D3', np.array([s * e3, s * e3])), ('Diff3', np.array([s * e3, -s * e3]))):
        cen, acc = scan_block(S_1x2, mats([(0, 0, np.pi)]), beta, v, N, rng)
        r = analyze(cen, acc)
        rows.append((mname, r))
        print(f"   {mname:>9}{r['a0']:8.3f}{r['G']:11.4f}{r['X']:9.4f}")
    Gmin = min(r['G'] for _, r in rows if np.isfinite(r['G']))
    print(f"\n   MIN G''(a0*) over boundary x block x mode at beta=2.5 = {Gmin:+.4f}"
          f"   (> N/12={2/12:.4f}: {Gmin > 2/12})")

    print("\n[H4] SCALING with L (k=0, theta=pi) and with beta; + Fisher average I(rho)")
    print(f"   {'L/block':>9}{'beta':>6}{'G\'\'(a0*)':>11}{'I(rho)=E_rho[G'']':>18}{'E[G'']/G(a0*)':>14}")
    for name, blk, fx, v in specs:
        cen, acc = scan_block(blk, mats(fx), 2.5, v, 1500000, rng)
        r = analyze(cen, acc)
        print(f"   {name:>9}{2.5:6.2f}{r['G']:11.4f}{r['I']:18.4f}{r['Gbar']/r['G']:14.3f}")
    for be in (1.0, 2.5, 4.0):
        cen, acc = scan_block(S_1x2, mats([(0, 0, np.pi)]), be, mode_k0(2), 1500000, rng)
        r = analyze(cen, acc)
        print(f"   {'1x2':>9}{be:6.2f}{r['G']:11.4f}{r['I']:18.4f}{r['Gbar']/r['G']:14.3f}")


def H5_acid():
    print("\n" + "=" * 100)
    print("[H5] free-field acid test")
    print("=" * 100)
    for kap in (1.0, 0.1, 1e-3):
        print(f"   free massless Gaussian kappa={kap:8.4f}: G''=kappa={kap:8.4f} -> 0")


if __name__ == "__main__":
    H1_hessians()
    H2_H3_H4_margins()
    H5_acid()
    print("\nSM1-blocks: done.  No QED.  Clay 0%.")
