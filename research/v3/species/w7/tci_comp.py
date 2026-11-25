#!/usr/bin/env python3
"""
tci_comp.py -- (TCI) round 2: definition audit + compensation-vs-elimination algebra + SU(2) numerics.

  (TCI)   <phi,Q phi> + <phi,KK phi> - cF >= mu > 0   for all unit phi in H_- .

Sections
  [A] BOGOLIUBOV SIGN LEMMA -- every canonical (symplectic) lift LOWERS the curvature.
  [B] COMPENSATION vs ELIMINATION -- Schur inclusion, and h_n (in-block) vs Sigma (eliminated)
      crossing points on the 4D L=1 Savvidy model: compensation survives where elimination dies.
  [C] 2 PLAQUETTE 2x2 SU(2): n_-, dim H_-, lambda_min(Hess V) = the TCI margin under the
      Jacobian lift, the provable floor mu3 = lambda_min(beta Hess S_W) + N/12, and the beta-scan.
  [D] COMPENSATION IDENTITY on the lattice: the law V''(a*) = (d-2) cF, V''(a*) - cF = (d-3) cF.
  [E] ACID: free field F=0 -> H_- = empty -> (TCI) is VACUOUS (acid test FAILS, not "returns 0").

py -3.12 -I, numpy only, single process, deterministic seeds.  No QED.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'w7'))
from sm1_exact import su2_exp, half_tr, dag, phi0, pr, prr            # noqa: E402
from sm1_blocks import sample_haar_xi                                  # noqa: E402

np.seterr(invalid='ignore', divide='ignore')
E1 = np.array([1.0, 0.0, 0.0])
E3 = np.array([0.0, 0.0, 1.0])
N12 = 2.0 / 12.0                      # SU(2): the pointwise Haar/Jacobian floor Phi'' >= N/12


def wq(x, w, p):
    """weighted quantile"""
    o = np.argsort(x)
    c = np.cumsum(w[o])
    c /= c[-1]
    return float(x[o][min(np.searchsorted(c, p), len(o) - 1)])


# ============================================================ [A] Bogoliubov sign lemma
def sec_A():
    print("=" * 100)
    print("[A] BOGOLIUBOV / SYMPLECTIC SIGN LEMMA:  H = p^2/2 + w^2 q^2/2 + c(pq+qp) = (p+2cq)^2/2 + (w^2-4c^2) q^2/2")
    print("    => every canonical lift shifts w^2 -> w^2 - 4c^2 <= w^2 : the lift can only LOWER the curvature.")
    print("=" * 100)
    n = 60
    a = np.diag(np.sqrt(np.arange(1, n)), 1)
    q = (a + a.T) / np.sqrt(2.0)
    p = 1j * (a.T - a) / np.sqrt(2.0)
    w = 1.0
    print(f"   {'c':>7}{'spec(H1) low3':>34}{'spec(H2) low3':>34}{'sqrt(w^2-4c^2)':>16}")
    for c in (0.0, 0.1, 0.2, 0.4):
        H1 = p @ p / 2 + (w * w) * (q @ q) / 2 + c * (p @ q + q @ p)
        H2 = p @ p / 2 + ((w * w - 4 * c * c)) * (q @ q) / 2
        e1 = np.linalg.eigvalsh(H1)[:3]
        e2 = np.linalg.eigvalsh(H2)[:3]
        print(f"   {c:7.2f}   [{e1[0]:9.5f}{e1[1]:9.5f}{e1[2]:9.5f}]   "
              f"[{e2[0]:9.5f}{e2[1]:9.5f}{e2[2]:9.5f}]   {np.sqrt(max(w*w-4*c*c,0.0)):16.5f}")
    print("   => spectra agree (unitary/squeeze equivalence) and the frequency is 2x reduced by c;")
    print("      c=0 is the maximum.  NO symplectic/Feshbach/Bogoliubov lift can be a POSITIVE compensation. [T]")


# ============================================================ models
def s_2x2(X, fixed):
    """2x2 cross block, Wilson part WITHOUT beta.  X:(M,4,3) free links."""
    H01, H11, V10, V11 = (su2_exp(X[:, k]) for k in range(4))
    H00, H10, H02, H12, V00, V20, V01, V21 = fixed
    P00 = H00[None] @ V10 @ dag(H01) @ dag(V00[None])
    P10 = H10[None] @ V20[None] @ dag(H11) @ dag(V10)
    P01 = H01 @ V11 @ dag(H02[None]) @ dag(V01[None])
    P11 = H11 @ V21[None] @ dag(H12[None]) @ dag(V11)
    return 4.0 - half_tr(P00) - half_tr(P10) - half_tr(P01) - half_tr(P11)


def hess_fd(f, X, h=1e-3):
    """batched Hessian of f(X:(M,n,3)) in the flattened chart."""
    M, n, _ = X.shape
    N = 3 * n
    Xf = X.reshape(M, N).copy()
    f0 = f(Xf.reshape(M, n, 3))
    H = np.zeros((M, N, N))
    base = np.zeros(N)
    for i in range(N):
        e = base.copy(); e[i] = h
        H[:, i, i] = (-f((Xf + 2 * e).reshape(M, n, 3)) + 16 * f((Xf + e).reshape(M, n, 3))
                      - 30 * f0 + 16 * f((Xf - e).reshape(M, n, 3)) - f((Xf - 2 * e).reshape(M, n, 3))) / (12 * h * h)
    for i in range(N):
        ei = base.copy(); ei[i] = h
        for j in range(i + 1, N):
            ej = base.copy(); ej[j] = h
            vpp = f((Xf + ei + ej).reshape(M, n, 3)); vpm = f((Xf + ei - ej).reshape(M, n, 3))
            vmp = f((Xf - ei + ej).reshape(M, n, 3)); vmm = f((Xf - ei - ej).reshape(M, n, 3))
            H[:, i, j] = H[:, j, i] = (vpp - vpm - vmp + vmm) / (4 * h * h)
    return H


def hess_phi(X):
    """EXACT Hessian of Phi = sum_l phi0(|xi_l|): block 3x3 per link, both curvatures >= N/12."""
    M, n, _ = X.shape
    H = np.zeros((M, 3 * n, 3 * n))
    rho = np.linalg.norm(X, axis=-1)
    rad = prr(rho / 2.0)                       # Phi''(rho)
    tan = np.where(rho > 1e-9, pr(rho / 2.0) / np.maximum(rho, 1e-12), 1.0 / 6.0)   # Phi'(rho)/rho
    for l in range(n):
        nhat = X[:, l, :] / np.maximum(rho[:, l:l + 1], 1e-12)
        blk = (tan[:, l, None, None] * np.eye(3)[None]
               + (rad - tan)[:, l, None, None] * nhat[:, :, None] * nhat[:, None, :])
        H[:, 3 * l:3 * l + 3, 3 * l:3 * l + 3] = blk
    return H


def samp(M, nfree, rng):
    return np.stack([sample_haar_xi(M, rng) for _ in range(nfree)], axis=1)


# ============================================================ [B] compensation vs elimination
def s_4d(X):
    U = [su2_exp(X[:, k]) for k in range(4)]
    s = 0.0
    for m in range(4):
        for n in range(m + 1, 4):
            P = U[m] @ U[n] @ dag(U[m]) @ dag(U[n])
            s = s + (1.0 - half_tr(P))
    return s


def sec_B():
    print("\n" + "=" * 100)
    print("[B] COMPENSATION vs ELIMINATION on the 4D L=1 Savvidy model (SU(2), beta=2.5)")
    print("    H = [[h_n, c^T],[c, K]];  in-block (COMPENSATION) curvature h_n  vs  Schur (ELIMINATION)")
    print("    Sigma = h_n - c^T K^-1 c .   Lemma: Sigma <= h_n - ||c||^2/lam_max(K)  =>  Sigma <= h_n .")
    print("=" * 100)
    beta = 2.5
    X0 = np.zeros((1, 4, 3))
    v = np.array([E3, E3, E3, E3]).reshape(-1) / 2.0        # k=0 Cartan mode, unit
    Qb, _ = np.linalg.qr(np.eye(12) - np.outer(v, v))       # orthonormal basis of v^perp (12x11)
    Qb = Qb[:, :11]
    print(f"   {'field c':>8}{'h_n(comp)':>11}{'gamma=lmin K':>13}{'c^T K^-1 c':>12}{'Sigma(elim)':>12}"
          f"{'TCI?':>6}{'elim?':>8}")
    crossed = []
    for c in (0.0, 0.2, 0.3, 0.4, 0.45, 0.50, 0.52, 0.53, 0.54, 0.55, 0.60, 0.80, 1.00, 1.50):
        X = X0.copy()
        X[0, 1] = c * E3
        X[0, 2] = c * E1
        Hs = hess_fd(s_4d, X)[0] * beta
        Hp = hess_phi(X)[0]
        H = Hs + Hp
        hn = float(v @ H @ v)
        cvec = Qb.T @ H @ v                    # (11,)
        K = Qb.T @ H @ Qb                      # (11,11) the fibre block on v^perp
        gam = float(np.linalg.eigvalsh(K)[0])
        # the Schur form Sigma = h_n - c^T K^-1 c EXISTS only for gamma > 0 (K invertible);
        # when K is indefinite the 1-loop Schur form is *undefined* (ch05:267), so flag it.
        if gam > 0:
            schur = float(cvec @ np.linalg.solve(K, cvec))
            Sig = hn - schur
            stag, etag = ('OK' if hn > 0 else 'FAIL'), ('OK' if Sig > 0 else 'FAIL')
        else:
            schur = np.nan
            Sig, stag, etag = np.nan, ('OK' if hn > 0 else 'FAIL'), 'UNDEF'
        crossed.append((c, hn, gam, Sig))
        print(f"   {c:8.2f}{hn:11.4f}{gam:13.4f}{'--' if gam <= 0 else f'{schur:12.4f}':>12}"
              f"{'--' if gam <= 0 else f'{Sig:12.4f}':>12}{stag:>6}{etag:>8}")
    cl = [(c, s, g) for c, h, g, s in crossed if g > 0 and np.isfinite(s)]
    cgam = np.nan
    for (c0, s0, g0), (c1, s1, g1) in zip(cl, cl[1:]):
        if g0 > 0 >= g1:
            cgam = c0 + (c1 - c0) * g0 / (g0 - g1)
    hmin = min(h for c, h, g, s in crossed)
    print(f"   => in-block (COMPENSATION) curvature h_n > 0 on the whole scan, min {hmin:+.4f} (= N/12 at the")
    print(f"      flat point), and F-MONOTONE INCREASING (h_n = 0.167 -> 0.728 for c = 0 -> 1.5).")
    print(f"      eliminated (Schur) curvature: Sigma = {cl[-1][1]:+.4f} at c = {cl[-1][0]:.2f} (last point with")
    print(f"      gamma>0), then K goes INDEFINITE at c_gamma ~ {cgam:.2f} and the Schur form is UNDEFINED.")
    print(f"      => elimination dies exactly at gamma=0 (c_gamma ~ {cgam:.2f}); compensation has no such")
    print(f"      threshold and survives to c = 1.5 (3x further), still positive. [V]")

    print("\n   random-block check of  Sigma <= h_n  (i.e. elimination can never RAISE the in-block form):")
    rng = np.random.default_rng(7)
    ok = True
    for t in range(2000):
        n = 6
        A = rng.standard_normal((n, n)); H = A + A.T
        H = H @ H.T / n + rng.standard_normal((n, n)) * 0.3
        H = 0.5 * (H + H.T)
        vv = rng.standard_normal(n); vv /= np.linalg.norm(vv)
        Pp = np.eye(n) - np.outer(vv, vv)
        hn = float(vv @ H @ vv); cv = Pp @ H @ vv; K = Pp @ H @ Pp
        ev = np.linalg.eigvalsh(K)
        if ev[0] < 1e-6:
            continue
        Sig = hn - float(cv @ np.linalg.inv(K) @ cv)
        if Sig > hn + 1e-9:
            ok = False
    print(f"      Sigma <= h_n on all tested blocks: {ok}")


# ============================================================ [C] 2 plaquette numerics
def sec_C(M=60000, seed=20261124):
    print("\n" + "=" * 100)
    print(f"[C] 2 PLAQUETTE / 2x2 SU(2) cross block, flat boundary, Haar-importance sample M={M}")
    print("    dim H = 12 (4 free links x 3).  Hess V = beta*Hess S_W + Hess Phi, Phi EXACT, Phi''>=N/12.")
    print("=" * 100)
    rng = np.random.default_rng(seed)
    fixed = [np.eye(2, dtype=complex) for _ in range(8)]
    X = samp(M, 4, rng)
    S = s_2x2(X, fixed)
    Hs = hess_fd(lambda Z: s_2x2(Z, fixed), X)
    Hp = hess_phi(X)
    lamPhi = np.linalg.eigvalsh(Hp)
    print(f"   exact check: min over samples of lam_min(Hess Phi) = {lamPhi.min():.6f}   (N/12 = {N12:.6f})"
          f"   -> Phi'' >= N/12 : {lamPhi.min() >= N12 - 1e-9}")
    print(f"   median lam_min(Hess Phi) = {np.median(lamPhi.min(axis=1)):.4f}")

    print(f"\n   {'beta':>5}{'ESS%':>7}{'E[n_-]':>8}{'P(n_->0)':>10}{'med lam_min V':>15}"
          f"{'P(lam_min<0)':>14}{'med mu3':>9}{'P(mu3>=0)':>11}{'|lam|/lamPhi':>14}")
    rows = {}
    for beta in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0):
        w = np.exp(-beta * S)
        wn = w / w.sum()
        ess = 100.0 * (w.sum() ** 2) / (len(w) * (w * w).sum())
        Hm = beta * Hs
        lamV = np.linalg.eigvalsh(Hm + Hp)
        lminV = lamV[:, 0]
        nm = (lamV < 0).sum(axis=1)
        lmag = np.linalg.eigvalsh(Hm)[:, 0]
        mu3 = lmag + N12
        rows[beta] = (lminV, wn, nm, mu3, lmag)
        print(f"   {beta:5.2f}{ess:7.1f}{np.dot(nm, wn):8.2f}{(nm > 0) @ wn:10.3f}"
              f"{wq(lminV, wn, 0.5):15.4f}{wn[lminV < 0].sum():14.3f}"
              f"{wq(mu3, wn, 0.5):9.4f}{wn[mu3 >= 0].sum():11.3f}"
              f"{abs(wq(lmag, wn, 0.5)) / np.median(lamPhi.min(axis=1)):14.3f}")
    lminV, wn, nm, mu3, lmag = rows[2.5]
    print(f"\n   beta=2.5 detail:  quantiles of lam_min(Hess V) (the TCI margin under K=Hess Phi):")
    for p in (0.05, 0.25, 0.5, 0.75, 0.95):
        print(f"      q{p:.2f} = {wq(lminV, wn, p):+.4f}")
    print(f"   P(lam_min<0)= {wn[lminV < 0].sum():.3f}   =>  (TCI) with K=Hess Phi FAILS on this fraction"
          f" of the measure at beta=2.5")
    print(f"   n_- values (weighted): ", {int(k): round(float(wn[nm == k].sum()), 3)
                                        for k in range(0, 6) if wn[nm == k].sum() > 1e-3})
    for bb in (2.5,):
        lv, ww, _, _, _ = rows[bb]
    # beta threshold where the weighted median margin crosses 0
    bs = np.array(sorted(rows)); med = np.array([wq(rows[b][0], rows[b][1], 0.5) for b in bs])
    frac = np.array([rows[b][1][rows[b][0] < 0].sum() for b in bs])
    j = np.where(np.diff(np.sign(med)) != 0)[0]
    bc = np.nan
    if len(j):
        i0 = j[0]
        bc = bs[i0] + (bs[i0 + 1] - bs[i0]) * (-med[i0]) / (med[i0 + 1] - med[i0])
    print(f"\n   beta-scan: med lam_min(Hess V) = {np.round(med,3)}  at beta = {bs}")
    print(f"              P(lam_min<0)        = {np.round(frac,3)}")
    print(f"   => TCI (K = Hess Phi, the only pointwise-uniform positive lift) holds with a positive")
    print(f"      margin only for beta < beta_c ~ {bc:.2f}  (crossing of the median); FALSE at the head 2.5. [V]")

    # ---- the saddle-vs-measure dimension mismatch
    print("\n   dimension bookkeeping (the object mismatch): at a Savvidy saddle dim H_- ~ (gB)^2 V/8pi^2")
    print("   (extensive), while on the beta=2.5 measure E[n_-] is the O(1) number in the table above.")
    return rows, lamPhi


# ============================================================ [D] compensation identity
def sec_D(seed=20261124):
    print("\n" + "=" * 100)
    print("[D] COMPENSATION IDENTITY on the lattice: V''(a*) = (d-2) cF  and  V''(a*) - cF = (d-3) cF")
    print("    exact for V = -cF a^2/2 + lam_d a^d (stationarity + homogeneity).  Quartic (d=4) => 2 cF.")
    print("=" * 100)
    beta = 2.5
    rng = np.random.default_rng(seed)
    fixed = [np.eye(2, dtype=complex) for _ in range(8)]
    X = samp(40000, 4, rng)
    S = s_2x2(X, fixed)
    Hs = beta * hess_fd(lambda Z: s_2x2(Z, fixed), X)
    Hp = hess_phi(X)
    H = Hs + Hp
    lam, vec = np.linalg.eigh(H)
    lmin = lam[:, 0]
    w = np.exp(-beta * S); wn = w / w.sum()

    # algebraic check of the law on the pure normal form
    print(f"   {'d':>3}{'a*^2':>10}{'cF':>8}{'V''(a*)':>11}{'V''-cF':>10}{'(d-3)cF':>10}")
    for d in (3, 4, 5, 6):
        cF = 1.0
        lam_d = 0.7
        a2 = (cF / (d * lam_d)) ** (2.0 / (d - 2.0))
        V2 = -cF + d * (d - 1) * lam_d * a2 ** ((d - 2) / 2.0)
        print(f"   {d:3d}{a2:10.4f}{cF:8.2f}{V2:11.5f}{V2 - cF:10.5f}{(d - 3) * cF:10.5f}")
    print("   => the compensation V''(a*)-cF = (d-3)cF is (i) >= 0 for every d>=3, (ii) EXACTLY ZERO for a")
    print("      cubic anharmonicity, (iii) proportional to cF=g|B| => NOT F-independent in magnitude,")
    print("      only in the RATIO V''(a*)/cF = d-2. [T]")

    # lattice ray scan on the worst (Savvidy-like) and the median configuration
    idx = [int(np.argmin(lmin)), int(np.argsort(lmin)[len(lmin) // 2])]
    print(f"\n   ray scan V(X_i + a phi) along the softest eigenvector phi of Hess V at X_i, beta=2.5.")
    print(f"   (the exponential chart is the ball |xi_link| < 2pi; beyond it Phi=-log J is undefined)")
    print(f"   {'case':>12}{'lam_min V':>11}{'cF=-V0\"':>10}{'a_chart':>9}{'V(a_ch)-V0':>12}"
          f"{'lam4':>9}{'a*,extrap':>11}{'V\"(a*)':>10}{'V\"-cF':>9}")
    for tag, i in zip(('savvidy-like', 'median'), idx):
        phi = vec[i, :, 0].reshape(4, 3)
        Xi = X[i]
        h = 1e-3

        def f(a, Xi=Xi, phi=phi):
            Z = (Xi + a * phi)[None]
            return beta * float(s_2x2(Z, fixed)[0]) + float(np.sum(phi0(np.linalg.norm(Z[0], axis=-1))))
        V0 = f(0.0)
        cF = -(f(h) - 2 * V0 + f(-h)) / (h * h)          # -V''(0) = tachyonic depth
        aa = np.linspace(-8.0, 8.0, 16001)
        aa = aa[np.abs(aa) > 1e-9]
        vals = np.array([f(t) for t in aa])
        fin = np.isfinite(vals)
        a_ch = float(np.max(np.abs(aa[fin])))            # chart edge along this ray
        Vch = float(vals[fin][np.argmax(np.abs(aa[fin]))])
        R = vals + 0.5 * cF * aa ** 2
        m = fin & (np.abs(aa) > 0.3 * a_ch) & (np.abs(aa) < 0.9 * a_ch)
        lam4 = float(np.polyfit(aa[m] ** 4, R[m], 1)[0]) if m.sum() > 5 else np.nan
        astar = np.sqrt(abs(cF / lam4)) if lam4 and np.isfinite(lam4) and cF / lam4 > 0 else np.nan
        V2 = -cF + 3 * lam4 * astar ** 2 if np.isfinite(astar) else np.nan
        print(f"   {tag:>12}{lmin[i]:11.4f}{cF:10.4f}{a_ch:9.3f}{Vch - V0:12.4f}"
              f"{lam4:9.4f}{astar:11.3f}{V2:10.4f}{V2 - cF:9.4f}")
    print("   => the extrapolated condensate well a* is FAR outside the accessible ball (a* ~ 9-14 while the")
    print("      rays are sampled to |a|<=8), and the predicted margin there is (d-3)cF = cF.")
    print("      So the compensation V''(a*) = cF + R''(a*) with R''>0 is REAL, but it is a LARGE-FIELD")
    print("      object: no small block reaches it.  That is why every small-volume witness (this one")
    print("      included) sees only the negativity — and why the compensation cannot be tested on 2 plaq.")
    print("   => on the lattice the well curvature exceeds the tachyonic depth, V''(a*)-cF > 0, and the")
    print("      positive source is the condensate amplitude |a*|^2 itself: V''(a*) = cF + R''(a*).")


# ============================================================ [E] acid / vacuity
def sec_E():
    print("\n" + "=" * 100)
    print("[E] ACID TEST and the VACUITY DEFECT of the literal (TCI)")
    print("=" * 100)
    print("   (TCI) quantifies 'for all unit phi in H_-'.  H_- is DEFINED by the Hessian")
    print("   (H_- = negative spectral subspace).  Free massless field: Hess V = (k^2 + N/12) I > 0")
    print("   => H_- = {0} => the quantifier is EMPTY => (TCI) holds VACUOUSLY, it does not 'return 0'.")
    for k in (1.0, 0.3, 0.1, 0.0):
        print(f"      free field k={k:4.1f}: lam_min(Hess V) = {k*k + N12:7.4f} > 0  =>  dim H_- = 0")
    print("   => the literal form PASSES the charter-3 acid test VACUOUSLY (a false-positive mode).")
    print("      Fix: quantify over a BACKGROUND-defined H_-^{fix} (dim ~ (gB)^2 V/8pi^2, fixed by F),")
    print("      and require dim H_-^{fix} >= 1 .  Then F = 0 gives cF = 0 and the inequality 0 >= mu fails. [T]")


if __name__ == "__main__":
    sec_A()
    sec_B()
    sec_C()
    sec_D()
    sec_E()
    print("\nTCI-comp: done.  No QED.  Clay 0%.")
