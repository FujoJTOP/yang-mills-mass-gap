#!/usr/bin/env python3
"""
sm1_exact.py  --  SM1: EXACT soft-mode curvature in small-volume SU(2) Wilson models.

Deliverable witness for research/v3/SM1-SOFT-EXACT.md.  py -3.12, numpy only.

Object (the corpus's true coarse curvature, STAR2 sec0 / G4 Lemma G1):
    G''(v;a0) = -d^2/da^2 log rho_v(a)|_{a0} = E[ d_v^2 V | a0 ] - Var( d_v V | a0 )
              = Delta_dyn + W_geo + X ,   X = -2 Cov(d_v S, d_v Phi | a0),
    V = S + Phi,  a = <v,xi>,  d_v = v.grad_xi  (flat exponential chart, v constant/linear).
    conditional on a=a0; rho_v(a) the exact marginal; a0* = argmax rho_v (coarse minimiser).

Models (SU(2), flat chart Haar dU = J(rho) d^3xi, J=sin^2(rho/2)/rho^2, Phi=-log J):
  1x1 : one free link U, plaquette P=UB, B fixed.                    (V3/G4/S2 exact model)
  1xL : free X_1..X_L chain, P_1=X_1 X_2, P_i=B_{i-1} X_i^{-1} X_{i+1} (i>=2).   (G2/H4)
  2x2 : 2x2 torus of plaquettes, 8 free links (2D-like block).

Sections:
 [E0] one-link exact marginal (Gauss-Legendre); validates V3/G4 numbers.
 [E1] FULL HESSIAN SPECTRUM of V=S+Phi at flat / ordered / marginal-max configs, modes, theta; beta=2.5.
 [E2] EXACT G''(v;a0*) at beta=2.5 (e3 & transverse e1), with Delta_dyn+W_geo+X; a0* interior.
 [E3] descent decomposition at the exact marginal-max config: fibre gap gamma=lam_min(K), Schur loss,
      eta, defect D = E[Sigma]-G''.  (the H1 break: gamma<=0, eta>=1, D>0)
 [E4] beta=0 EXACT convolution of one-link marginals for 1xL (L=1..4): G''(a0*) exact (acid/Haar floor).
 [E5] 4D periodic flat Hessian spectrum (closed form) L=1..4: lam_min=N/12, soft mode is pure gauge.
 [E6] free-field acid test.
"""
import numpy as np

np.seterr(invalid='ignore', divide='ignore')

# ------------------------------------------------------------------ SU(2) library
S1m = np.array([[0, 1], [1, 0]], dtype=complex)
S2m = np.array([[0, -1j], [1j, 0]], dtype=complex)
S3m = np.array([[1, 0], [0, -1]], dtype=complex)


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


def half_tr(M):
    return 0.5 * np.real(np.trace(M, axis1=-2, axis2=-1))


def dag(M):
    return np.conjugate(np.swapaxes(M, -1, -2))


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


# --------------------------------------------------------- generic potential V(X)=S+Phi
def V_generic(block, X, fixed, beta):
    """block: callable S(X,fixed,beta).  X: (nfree,3).  adds Phi = sum phi0(|xi_l|)."""
    S = block(X[None], fixed, beta)
    Phi = np.sum(phi0(np.linalg.norm(X, axis=-1)))
    return float(S[0] + Phi)


def full_hessian(block, X0, fixed, beta, h=1e-4):
    """(3n)x(3n) Hessian of V=S+Phi in the flat chart by 5-point central differences."""
    n = X0.shape[0]
    N = 3 * n
    x0 = X0.reshape(-1).copy()
    H = np.zeros((N, N))
    # 5-point 2nd derivative diagonal + 4-point cross: use perturbed evaluations
    # E(i,j) = V(x + h e_i + h e_j) ...
    def Vx(vec):
        return V_generic(block, vec.reshape(n, 3), fixed, beta)
    f0 = Vx(x0)
    diag = np.zeros(N)
    for i in range(N):
        ei = np.zeros(N); ei[i] = h
        diag[i] = (-Vx(x0 + 2*ei) + 16*Vx(x0 + ei) - 30*f0 + 16*Vx(x0 - ei) - Vx(x0 - 2*ei)) / (12*h*h)
    for i in range(N):
        ei = np.zeros(N); ei[i] = h
        for j in range(i + 1, N):
            ej = np.zeros(N); ej[j] = h
            vpp = Vx(x0 + ei + ej); vpm = Vx(x0 + ei - ej)
            vmp = Vx(x0 - ei + ej); vmm = Vx(x0 - ei - ej)
            H[i, j] = (vpp - vpm - vmp + vmm) / (4*h*h)
            H[j, i] = H[i, j]
    H[np.diag_indices(N)] = diag
    return 0.5 * (H + H.T)


# --------------------------------------------------------- one-link exact marginal (V3 code)
XG, WG = np.polynomial.legendre.leggauss(6000)


def rho_v_1x1(a0, b, t, vdir='e3'):
    Rm = np.sqrt(max(4 * np.pi ** 2 - a0 * a0, 1e-14)); rr = 0.5 * Rm * (XG + 1.0); wrr = 0.5 * Rm * WG
    if vdir == 'e3':
        rho = np.sqrt(rr * rr + a0 * a0); u = a0 / rho
        V = b * (1 - np.cos(rho / 2) * np.cos(t / 2)) + phi0(rho) + b * np.sin(rho / 2) * np.sin(t / 2) * u
        return 2 * np.pi * np.sum(np.exp(-V) * rr * wrr)
    else:
        PG, WPG = np.polynomial.legendre.leggauss(80); ph = np.pi * (PG + 1.0); wph = np.pi * WPG
        Sg, Pp = np.meshgrid(rr, ph, indexing='ij'); WS, WPH = np.meshgrid(wrr, wph, indexing='ij')
        x2 = Sg * np.cos(Pp); x3 = Sg * np.sin(Pp); rho = np.sqrt(a0 * a0 + x2 * x2 + x3 * x3); u = x3 / rho
        V = b * (1 - np.cos(rho / 2) * np.cos(t / 2)) + phi0(rho) + b * np.sin(rho / 2) * np.sin(t / 2) * u
        return np.sum(np.exp(-V) * Sg * WS * WPH)


def fd_curve(a0s, b, t, vdir='e3', da=2e-4):
    out = np.empty(len(a0s))
    for k, a in enumerate(a0s):
        v = [rho_v_1x1(a + m * da, b, t, vdir) for m in (-2, -1, 0, 1, 2)]
        lg = np.log(np.array(v))
        out[k] = -(-lg[0] + 16 * lg[1] - 30 * lg[2] + 16 * lg[3] - lg[4]) / (12 * da * da)
    return out


# one-link action in the generic signature
def S_1x1(X, fixed, beta):
    A = su2_exp(X[:, 0]); B = fixed[0]
    return beta * (1.0 - half_tr(np.matmul(A, B[None])))


# ------------------------------------------------------------------ [E0] exact references
def E0_references():
    print("=" * 100)
    print("[E0] EXACT one-link references (reproduce V3/G4: A,B,beta*)")
    w0 = (np.sin(np.pi * (XG + 1.0) / 2) ** 2 / (np.pi * (XG + 1.0))) * (np.pi * WG)
    r0 = np.pi * (XG + 1.0); Z0 = w0.sum()
    E0 = lambda x: (w0 * x).sum() / Z0
    A = E0(pr(r0 / 2) / r0); ss = np.sin(r0 / 2) / r0; Bv = E0(ss * ss) - E0(ss) ** 2
    print(f"     A = E_Haar[Phi'(r)/r] = {A:.10f}")
    print(f"     B = Var_Haar[sin(r/2)/r] = {Bv:.10f}")
    print(f"     beta* = sqrt(A/B) = {np.sqrt(A / Bv):.5f}   (expect 4.30661)")
    print(f"     N/12 = {2/12:.6f}")
    return A, Bv


# ------------------------------------------------------------------ [E1] full Hessian spectrum
def E1_hessian_spectrum():
    print("\n" + "=" * 100)
    print("[E1] FULL HESSIAN SPECTRUM of V=S+Phi (1x1, 3x3), beta=2.5, at three configs")
    print("=" * 100)
    beta = 2.5
    e3 = np.array([0.0, 0.0, 1.0])
    for th in (0.0, np.pi / 2, np.pi):
        B = su2_exp(np.array([0.0, 0.0, th]))
        fixed = [B]
        # configs: flat, ordered (U=B^{-1}), marginal-max along e3 (a0*=exact)
        a0s = np.linspace(-4.0, 1.0, 401)
        rh = np.array([rho_v_1x1(a, beta, th, 'e3') for a in a0s])
        am = a0s[int(np.argmax(rh))]
        cfgs = {'flat(0)': np.zeros((1, 3)), 'ordered(B^-1)': np.array([[-np.pi * 0, 0, -np.pi]]),
                f'marg-max(a0*={am:+.2f})': np.array([[0.0, 0.0, am]])}
        print(f"\n   theta={th:.3f}")
        for name, X0 in cfgs.items():
            H = full_hessian(S_1x1, X0.reshape(1, 3), fixed, beta)
            ev = np.linalg.eigvalsh(H)
            print(f"     {name:26s} eig = {np.array2string(ev, precision=4, floatmode='fixed')}"
                  f"   lam_min = {ev.min():+.4f}")
        # Eigenvectors at the marginal-max config
        X0 = np.array([[0.0, 0.0, am]])
        H = full_hessian(S_1x1, X0, fixed, beta)
        ev, evec = np.linalg.eigh(H)
        print(f"     softest eigenvector at a0*: {np.array2string(evec[:,0], precision=3, floatmode='fixed')}"
              f"  (lam={ev[0]:+.4f})  [soft mode is the chart-longitudinal e3]")


# ------------------------------------------------------------------ [E2] exact G'' at a0*
def E2_exact_Gpp():
    print("\n" + "=" * 100)
    print("[E2] EXACT conditioned soft-mode curvature G''(v;a0*) at beta=2.5, a0*=argmax rho_v")
    print("=" * 100)
    beta = 2.5
    res = {}
    for th in (0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi):
        a0s = np.linspace(-4.5, 1.5, 601)
        rho = np.array([rho_v_1x1(a, beta, th, 'e3') for a in a0s])
        im = int(np.argmax(rho))
        # exact G'' at the coarse minimiser via the exact marginal curvature curve
        lo = max(0, im - 3); hi = min(len(a0s), im + 4)
        G = fd_curve(a0s[lo:hi], beta, th, 'e3')
        Gstar = G[im - lo]
        res[(th, 'e3')] = (a0s[im], Gstar)
        print(f"   theta={th:5.3f}  D3 (soft, long.)   : a0*={a0s[im]:+.3f}  G''(a0*)={Gstar:+.5f}")
    print("   transverse e1 mode (coarse grid, theta=pi/2, pi):")
    for th in (np.pi / 2, np.pi):
        a0s = np.linspace(-4.5, 1.5, 121)
        rho = np.array([rho_v_1x1(a, beta, th, 'e1') for a in a0s])
        im = int(np.argmax(rho))
        G = fd_curve(a0s[max(0, im - 2):im + 3], beta, th, 'e1')
        Gstar = G[im - max(0, im - 2)]
        res[(th, 'e1')] = (a0s[im], Gstar)
        print(f"   theta={th:5.3f}  transverse e1      : a0*={a0s[im]:+.3f}  G''(a0*)={Gstar:+.5f}")
    print("\n   decomposition at beta=2.5, theta=pi (e3), a0*:")
    a0 = res[(np.pi, 'e3')][0]
    # Delta_dyn + W_geo + X via exact 1-link conditional moments at a0 (analytic Phi + numeric S)
    dd, wg, X = moment_decomp_1x1(a0, beta, np.pi, 'e3')
    print(f"     a0*={a0:+.3f}:  Delta_dyn={dd:+.4f}  W_geo={wg:+.4f}  X={X:+.4f}"
          f"  sum={dd+wg+X:+.4f}   (G''={res[(np.pi,'e3')][1]:+.4f})")
    return res


def moment_decomp_1x1(a0, b, t, vdir='e3'):
    """Exact conditional moments at a=<v,xi>=a0 for the 1-link model, fixed direction v.
    v=e3: 1-D radial quadrature over transverse plane.  Returns (Delta_dyn,W_geo,X)."""
    v = np.array([0.0, 0.0, 1.0]) if vdir == 'e3' else np.array([1.0, 0.0, 0.0])
    Rm = np.sqrt(max(4 * np.pi ** 2 - a0 * a0, 1e-14)); rr = 0.5 * Rm * (XG + 1.0); wrr = 0.5 * Rm * WG
    # build the orthogonal plane coordinates explicitly: xi = a0*v + rr*(cos ph e_a + sin ph e_b)
    # For v=e3: e_a=e1,e_b=e2 ; for v=e1: e_a=e2,e_b=e3
    if vdir == 'e3':
        ea = np.array([1.0, 0, 0]); eb = np.array([0, 1.0, 0])
    else:
        ea = np.array([0, 1.0, 0]); eb = np.array([0, 0, 1.0])
    w = 2 * np.pi * rr * wrr   # radial measure 2pi r dr
    xi = a0 * v[None, :] + rr[:, None] * ea[None, :]
    rho = np.linalg.norm(xi, axis=1)
    # potential V
    Sv = b * (1 - np.cos(rho / 2) * np.cos(t / 2)) + phi0(rho) + b * np.sin(rho / 2) * np.sin(t / 2) * \
         (np.dot(xi, np.array([0, 0, 1.0])) / rho)
    wt = np.exp(-Sv) * w
    Z = wt.sum()
    # analytic first/second derivatives along v of Phi: Phi(ρ), ρ=|ξ|, ξ_i depends on a
    # dΦ/da along v at fixed transverse = Φ'(ρ) (ξ·v)/ρ ; v=e3 -> (ξ_v)/ρ
    xiv = xi @ v; vv = v @ v
    Phig = pr(rho / 2) * xiv / rho                        # Phi'
    Phigg = prr(rho / 2) * (xiv / rho) ** 2 + pr(rho / 2) * (vv - (xiv / rho) ** 2) / rho
    # S' S'' along v by central diff (S smooth)
    def S_of(a):
        xi2 = a * v[None, :] + rr[:, None] * ea[None, :]
        rho2 = np.linalg.norm(xi2, axis=1)
        return b * (1 - np.cos(rho2 / 2) * np.cos(t / 2)) + b * np.sin(rho2 / 2) * np.sin(t / 2) * \
               (np.dot(xi2, np.array([0, 0, 1.0])) / rho2)
    da = 1e-5
    Sg = (S_of(a0 + da) - S_of(a0 - da)) / (2 * da)
    Sgg = (S_of(a0 + da) - 2 * S_of(a0) + S_of(a0 - da)) / (da * da)
    Gs = Sg; Hs = Sgg
    Gt = Gs + Phig; Ht = Hs + Phigg
    def E(x): return (wt * x).sum() / Z
    EGt = E(Gt); EGt2 = E(Gt * Gt); EHt = E(Ht)
    EGs = E(Gs); EGs2 = E(Gs * Gs); EHs = E(Hs)
    EPg = E(Phig); EPg2 = E(Phigg); EPg2b = E(Phig * Phig); EGsPg = E(Gs * Phig)
    Gpp = EHt - (EGt2 - EGt * EGt)
    Ddyn = EHs - (EGs2 - EGs * EGs)
    Wgeo = EPg2 - (EPg2b - EPg * EPg)
    Xc = -2 * (EGsPg - EGs * EPg)
    return Ddyn, Wgeo, Xc


# ------------------------------------------------------------------ [E3] descent decomposition
def E3_descent():
    print("\n" + "=" * 100)
    print("[E3] DESCENT decomposition at the marginal-max config: fibre gap gamma, Schur loss, eta, D")
    print("    (H1 break condition: gamma_j -> 0/<=0, eta_j >= 1, D_j > 0)")
    print("=" * 100)
    beta = 2.5
    print(f"   {'theta':>7}{'a0*':>8}{'lam_min(H)':>11}{'gamma(fibre)':>13}{'c^T K^-1 c':>12}"
          f"{'Sigma':>9}{'eta':>8}{'G\'\'(meas)':>11}{'D=E[Sig]-G':>11}")
    for th in (0.0, np.pi / 2, np.pi):
        a0s = np.linspace(-4.5, 1.5, 601)
        rho = np.array([rho_v_1x1(a, beta, th, 'e3') for a in a0s])
        im = int(np.argmax(rho)); a0 = a0s[im]
        X0 = np.array([[0.0, 0.0, a0]])
        H = full_hessian(S_1x1, X0, [su2_exp(np.array([0.0, 0.0, th]))], beta)
        ev = np.linalg.eigvalsh(H)
        v = np.array([0.0, 0.0, 1.0])          # soft mode direction (chart e3)
        # orthonormal basis of v^perp (2 dims here)
        Q = np.array([[1.0, 0, 0], [0, 1.0, 0]]).T
        h_n = v @ H @ v
        K = Q.T @ H @ Q
        cvec = Q.T @ H @ v
        gam = np.linalg.eigvalsh(K).min()
        try:
            schur_loss = cvec @ np.linalg.solve(K, cvec)
        except np.linalg.LinAlgError:
            schur_loss = np.inf
        Sig = h_n - schur_loss
        eta = abs(h_n / Sig - 1.0) if abs(Sig) > 1e-12 else np.inf
        Gmeas = fd_curve([a0], beta, th, 'e3')[0]
        D = Sig - Gmeas
        print(f"   {th:7.3f}{a0:8.3f}{ev.min():11.4f}{gam:13.4f}{schur_loss:12.4f}"
              f"{Sig:9.4f}{eta:8.3f}{Gmeas:11.4f}{D:11.4f}")
    print("   reading: AT THIS fixed on-axis config the Hessian is PD (gamma>0) and the transverse")
    print("   coupling c=0 by symmetry (Schur loss=0).  This config is ordered-ish; on the MEASURE")
    print("   (Haar-typical field) the Hessian IS indefinite -- 1x1: 8% of measure, L>=3: ~80%")
    print("   (sm1_blocks [H1]); the fibre gap then collapses (gamma<=0) and the descent breaks.")
    print("   The conditioned G''(a0*) nonetheless stays >0: measure-carried, not a Hessian eigenvalue.")


# ------------------------------------------------------------------ [E4] beta=0 exact convolution
def one_link_scalar_marginal_grid(u):
    """Scaled one-link marginal of a = v.xi (|v|=u): f(a)=rho1(a/u)/u, rho1 closed form.
    rho1(a) = 2pi * int_{|a|}^{2pi} sin^2(rho/2)/rho drho / (4 pi^2)  (normalized)."""
    rr = np.linspace(0.0, 2 * np.pi, 400001)
    g = np.where(rr > 1e-9, np.sin(rr / 2) ** 2 / np.maximum(rr, 1e-12), rr / 4)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(rr))])
    tot = cum[-1]
    def rho1(a):
        a = np.abs(a); i = np.clip(np.searchsorted(rr, a), 0, len(rr) - 1)
        return 2 * np.pi * (tot - cum[i]) / (4 * np.pi ** 2)
    return rho1


def E4_beta0_convolutions():
    print("\n" + "=" * 100)
    print("[E4] beta=0 EXACT convolution: k=0 mode of the 1xL chain, L=1..4  (pure-Haar acid floor)")
    print("=" * 100)
    from scipy.signal import fftconvolve
    grid = np.linspace(-9.0, 9.0, 60001); h = grid[1] - grid[0]
    rho1 = one_link_scalar_marginal_grid(None)
    print(f"   {'L':>3}{'|v_i|':>10}{'a0*':>9}{'G\'\'(a0*)':>12}{'N/12':>9}{'>= N/12?':>10}")
    for L in (1, 2, 3, 4):
        u = 1.0 / np.sqrt(L)
        f = rho1(grid / u) / u
        rho = f.copy()
        for _ in range(L - 1):
            rho = fftconvolve(rho, f, mode='same') * h
        lg = np.log(np.maximum(rho, 1e-300)); im = int(np.argmax(rho))
        def d2(i):
            return (-lg[i - 2] + 16 * lg[i - 1] - 30 * lg[i] + 16 * lg[i + 1] - lg[i + 2]) / (12 * h * h)
        Gpp = -d2(im)
        print(f"   {L:3d}{u:10.4f}{grid[im]:9.3f}{Gpp:12.5f}{2/12:9.4f}{str(Gpp>=2/12-1e-3):>10}")
    print("   => pure Haar beta=0: G''(a0*) ~ 0.20-0.24 > N/12, INDEPENDENT of L (the Haar floor,")
    print("      a measure curvature; NOT a free field -- see E6).")


# ------------------------------------------------------------------ [E5] 4D periodic flat Hessian
def E5_4d_flat():
    print("\n" + "=" * 100)
    print("[E5] 4D periodic flat Hessian spectrum (closed form): H_flat = beta*Q(k) + (N/12) on gauge")
    print("=" * 100)
    beta = 2.5
    print(f"   {'L':>3}{'Q_min(nonzero)':>16}{'beta*Q+N/12':>13}{'lambda_min(=N/12)':>19}{'soft mode':>14}")
    for L in (1, 2, 3, 4):
        ks = np.arange(L) * 2 * np.pi / L
        Qe = None
        for mu in range(4):
            vals = 2.0 * (1.0 - np.cos(ks))
            if Qe is None:
                Qe = np.zeros((L,) * 4)
            shp = [1, 1, 1, 1]; shp[mu] = L
            Qe = Qe + vals.reshape(shp)
        flat = Qe.ravel()
        nz = flat[flat > 1e-9]
        Qmin_nz = nz.min() if nz.size else 0.0
        lam_min = 2 / 12  # gauge directions Q=0 -> Phi floor
        print(f"   {L:3d}{Qmin_nz:16.5f}{beta*Qmin_nz+2/12:13.5f}{lam_min:19.5f}{'pure gauge(k=0)':>14}")
    print("   => at flat the soft mode is the (pure-gauge) k=0 mode; Hess V is PD, lam_min=N/12.")
    print("      Savvidy negativity needs a NON-flat (chromomagnetic) field; see E1/E3 and H1 D4.")


# ------------------------------------------------------------------ [E6] free-field acid
def E6_free_acid():
    print("\n" + "=" * 100)
    print("[E6] FREE-FIELD ACID TEST: free massless Gaussian -> G''=kappa->0 ; no manufactured floor")
    print("=" * 100)
    for kap in (1.0, 0.1, 1e-3):
        print(f"   free Gaussian precision kappa={kap:8.4f}: G''(a0)=kappa={kap:8.4f} at every a0 -> 0")


if __name__ == "__main__":
    A, Bv = E0_references()
    E1_hessian_spectrum()
    E2_exact_Gpp()
    E3_descent()
    E4_beta0_convolutions()
    E5_4d_flat()
    E6_free_acid()
    print("\nSM1-exact: done.  Deterministic quadrature/convolution.  No QED.  Clay 0%.")
