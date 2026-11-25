"""mc_ricci.py -- N1 (metric route).  Can the choice of Riemannian metric g flip the sign
of the BD/Polchinski criterion (star)?

All geometry is EXACT (Koszul, no finite differences).

FRAME.  e_a = lam_a = i T_a, tr(T_a T_b) = delta_ab/2, so g(X,Y) = -2 tr(XY) makes {lam_a}
orthonormal.  delta_ab = "the" Ad-invariant (bi-invariant) form.

LEFT-INVARIANT METRIC g_G(X,Y) = sum_ab G_ab x^a y^b, G symmetric > 0.
  * Ad-invariance <=> G = c delta  (su(N) simple: the invariant inner product is unique up
    to scale).  So all shape freedom of this route is NON-Ad-invariant G.
  * G -> cG: the Levi-Civita connection is unchanged, hence Hess_{cG} = Hess_G and
    Ric_{cG} = Ric_G as (0,2) tensors, while  (star) "Hess V + Ric/2 >= lam_dot G"
    becomes lam_dot(cG) = lam_dot(G)/c.  The final LSI carries |grad F|^2_{C_0}, C_0 = G^{-1},
    which scales the other way: the LSI is EXACTLY invariant.  => the SCALE of g is not free
    data; only the SHAPE (projective class) is.
  * FRAME INVARIANCE (used to cut the search from n(n+1)/2 to n params): an orthogonal frame
    change O in O(n) fixes delta and carries G -> O G O^T, Ric_G -> O Ric_G O^T,
    Hess_G -> O Hess_G O^T.  Hence rho(G,W) below is invariant, and WLOG G is DIAGONAL.
    The shape space is therefore effectively n eigenvalues (mod overall scale).  Verified
    numerically at the end (random non-diagonal G never beats the diagonal optimum).

KOSZUL (left-invariant, all derivative terms vanish):
  f^c_ab = structure constants,  f_abc = G_cd f^d_ab
  Gamma_{c,ab} = (1/2)(f_abc - f_acb - f_bca),  Gamma^c_ab = G^cd Gamma_{d,ab}
  Ric_ab = Gamma^c_cd Gamma^d_ab - Gamma^c_bd Gamma^d_ac
  Hess_G phi(e_a,e_b) = e_a e_b phi - Gamma^c_ab e_c phi

SINGLE PLAQUETTE, phi(U) = Re tr(U W), W := U2 U3 U4 U1 (the plaquette holonomy), evaluated
at U1 = 1 in the left frame:
  e_a e_b phi = Re tr(e_a e_b W),  e_c phi = Re tr(e_c W)
  Hess S_p = -(beta/N)[ D^2 phi - Gamma^c_ab (e_c phi) ]
One link sits in 6 plaquettes in d=4, so
  M_G(W) = 6 Hess S_p(W) + Ric_G/2,   rho(G,W) = lam_min(G^{-1/2} M_G G^{-1/2}).
rho > 0 <=> M_G is a positive definite (0,2)-tensor (congruence invariance => sign of rho is
a property of the tensor, independent of how G normalises it).

Q1  sup_G lam_min(Ric_G)  =?  N/4     (is the bi-invariant metric already optimal?)
Q2  sup_G inf_W rho(G,W)  : can a metric shape make lambda_dot_0 > 0 ?
"""
import numpy as np
import json
import scipy.linalg
from scipy.optimize import minimize

np.set_printoptions(precision=6, suppress=True, linewidth=160)
rng = np.random.default_rng(20261125)
BETA = 2.5


def su_generators(N):
    Ts = []
    for j in range(N):
        for k in range(j + 1, N):
            T = np.zeros((N, N), dtype=complex); T[j, k] = T[k, j] = 0.5; Ts.append(T)
    for j in range(N):
        for k in range(j + 1, N):
            T = np.zeros((N, N), dtype=complex); T[j, k] = -0.5j; T[k, j] = 0.5j; Ts.append(T)
    for l in range(1, N):
        T = np.zeros((N, N), dtype=complex)
        d = 1.0 / np.sqrt(2.0 * l * (l + 1))
        for j in range(l):
            T[j, j] = d
        T[l, l] = -l * d
        Ts.append(T)
    return Ts


def frame(N):
    lam = [1j * T for T in su_generators(N)]
    n = len(lam)
    f = np.zeros((n, n, n))
    for a in range(n):
        for b in range(n):
            C = lam[a] @ lam[b] - lam[b] @ lam[a]
            f[:, a, b] = [-2.0 * np.trace(C @ lam[c]).real for c in range(n)]
    return lam, f


def gamma_of(G, f):
    """Gamma^c_ab with index layout [c,a,b] = Gamma_{c,ab} raised.
    T[i,j,k] = f_{ijk}.  Gamma_{c,ab} = (1/2)(f_abc - f_acb - f_bca), so as an array indexed
    [c,a,b] the three terms are  T.transpose(1,2,0) , T.transpose(1,0,2) , T.transpose(2,0,1)
    (check: T.transpose(1,2,0)[c,a,b] = T[a,b,c] = f_abc, etc.)."""
    Gi = np.linalg.inv(G)
    T = np.einsum('cd,dab->abc', G, f)          # T[a,b,c] = f_abc = G_cd f^d_ab
    Glow = 0.5 * (np.transpose(T, (1, 2, 0)) - np.transpose(T, (1, 0, 2))
                  - np.transpose(T, (2, 0, 1)))
    return np.einsum('cd,dab->cab', Gi, Glow)


def ricci_of(G, f):
    Gm = gamma_of(G, f)
    return (np.einsum('ccd,dab->ab', Gm, Gm) - np.einsum('cbd,dac->ab', Gm, Gm))


def phi_data(W, lam):
    """(D^2 phi_ab, e_c phi) at U = 1 for phi(U) = Re tr(U W).  G-independent."""
    lw = np.array([np.trace(l @ W).real for l in lam])
    D2 = np.array([[np.trace(lam[a] @ lam[b] @ W).real for b in range(len(lam))]
                   for a in range(len(lam))])
    return D2, lw


def torus_elem(N, th):
    d = np.concatenate([np.exp(1j * np.asarray(th, float)), [np.exp(-1j * np.sum(th))]])
    return np.diag(d)


def rho_of(G, data, f, N, beta=BETA, npla=6.0):
    """lam_min(G^{-1/2} M G^{-1/2}) via the generalized symmetric problem M v = mu G v."""
    D2, lw = data
    H = D2 - np.einsum('cab,c->ab', gamma_of(G, f), lw)
    M = -npla * (beta / N) * H + 0.5 * ricci_of(G, f)
    M = 0.5 * (M + M.T)
    return scipy.linalg.eigh(M, G, eigvals_only=True)[0]


def diag_from_logx(x):
    """det-normalised positive diagonal matrix from n-1 free logs (clipped: the shape space is
    non-compact, and lambda_dot -> -inf at the boundary, so clipping is harmless for a max)."""
    n = len(x) + 1
    lx = np.concatenate([np.asarray(x, float), [-np.sum(x)]])
    return np.diag(np.exp(np.clip(lx, -12.0, 12.0)))


print("=" * 104)
print("STEP 1  sanity: Koszul Gamma(delta) = f/2 and Ric(delta) = (N/4) delta")
for N in (2, 3):
    lam, f = frame(N)
    n = len(lam)
    eg = np.abs(gamma_of(np.eye(n), f) - 0.5 * f).max()
    R = ricci_of(np.eye(n), f)
    print(f"  SU({N}) dim={n:2d}  max|Gamma(delta)-f/2|={eg:.2e}  "
          f"max|Ric(delta)-(N/4)delta|={np.abs(R - N/4*np.eye(n)).max():.2e}")

print()
print("=" * 104)
print("STEP 2  (T1) explicit SU(2) 3x3 toy: how Hess depends on G   [W = -1, beta=2.5]")
N = 2
lam2, f2 = frame(2)
Wm = -np.eye(2, dtype=complex)
d2 = phi_data(Wm, lam2)
for name, g in [("delta (bi-invariant)", [1, 1, 1]), ("diag(1,1,2)", [1, 1, 2]),
                ("diag(1,1,.5)", [1, 1, .5]), ("diag(1,4,4)", [1, 4, 4])]:
    G = np.diag(np.array(g, float))
    H = d2[0] - np.einsum('cab,c->ab', gamma_of(G, f2), d2[1])
    R = ricci_of(G, f2)
    print(f"  G={name:22s} 6*Hess S_p = {np.diag(-6*(BETA/2)*H)}   "
          f"Ric_G(diag) = {np.diag(R)}   eig(Ric) = {np.linalg.eigvalsh(R)}   "
          f"rho(G,-1) = {rho_of(G, d2, f2, 2):+.5f}")
print("  closed form check (SU(2), diagonal G): Ric_11 = [g1^2-(g2-g3)^2]/(2 g2 g3), cyclic")
for g in ([1, 1, 1], [1, 1, 2], [1, 4, 4]):
    g1, g2, g3 = g
    cf = [(g1**2 - (g2 - g3)**2) / (2 * g2 * g3), (g2**2 - (g3 - g1)**2) / (2 * g3 * g1),
          (g3**2 - (g1 - g2)**2) / (2 * g1 * g2)]
    print(f"     G=diag{g}: closed form {np.round(cf,6)}  code {np.round(np.diag(ricci_of(np.diag(np.array(g,float)),f2)),6)}")

print()
print("=" * 104)
print("STEP 3  (Q1) sup over SHAPES of lam_min(Ric_G)   [delta = I, det G = 1]")
res = {"beta": BETA}
for N in (2, 3):
    lam, f = frame(N)
    n = len(lam)

    def obj(x):
        return np.linalg.eigvalsh(ricci_of(diag_from_logx(x), f))[0]

    best, bx = -1e9, None
    for s in range(30):
        x0 = np.zeros(n - 1) if s == 0 else rng.normal(0, 1.5, n - 1)
        r = minimize(lambda x: -obj(x), x0, method='Nelder-Mead',
                     options=dict(maxiter=6000, fatol=1e-13, xatol=1e-10))
        if -r.fun > best:
            best, bx = -r.fun, r.x.copy()
    Gb = diag_from_logx(bx)
    # random NON-diagonal shapes, to confirm the diagonal WLOG
    rnd_best = -1e9
    for _ in range(4000):
        S = rng.normal(0, 1.2, (n, n)); S = (S + S.T) / 2
        w, V = np.linalg.eigh(S)
        G = V @ np.diag(np.exp(w - w.mean())) @ V.T
        rnd_best = max(rnd_best, np.linalg.eigvalsh(ricci_of(G, f))[0])
    print(f"  SU({N}): sup_G lam_min(Ric_G) = {best:.6f}    N/4 = {N/4:.6f}    "
          f"{'SAME (bi-invariant optimal)' if abs(best - N/4) < 2e-3 else 'DIFFERENT!'}")
    print(f"          eigen(Ric at optimum) = {np.linalg.eigvalsh(ricci_of(Gb, f))}")
    print(f"          best of 4000 random NON-diagonal G = {rnd_best:.6f}  (<= {best + 2e-3})")
    res[f"q1_SU{N}"] = dict(sup=float(best), biinv=float(N / 4), random_non_diag=float(rnd_best))

def mu_of(G, data, f, N, beta=BETA, npla=6.0):
    """lambda_min of the (0,2)-TENSOR M_G = 6 Hess_G S_p + Ric_G/2, evaluated as a matrix in
    the bi-invariant frame delta = I.  Scale invariant in G: M_{cG} = M_G.  Because
    lambda_min(G^{-1/2} M G^{-1/2}) and M have the same INERTIA (congruence), sign(mu) is
    exactly the sign of the criterion's operator: mu > 0 <=> M_G is positive definite."""
    D2, lw = data
    H = D2 - np.einsum('cab,c->ab', gamma_of(G, f), lw)
    M = -npla * (beta / N) * H + 0.5 * ricci_of(G, f)
    return np.linalg.eigvalsh(0.5 * (M + M.T))[0]


print()
print("=" * 104)
print("STEP 4  (Q2) inf_W sup_G mu(G,W): is there a configuration that NO metric shape can fix?")
for N in (2, 3):
    lam, f = frame(N)
    n = len(lam)
    if N == 2:
        Ws = [torus_elem(2, [t]) for t in np.linspace(0, 1, 201) * np.pi]
    else:
        gr = np.linspace(-1, 1, 17) * np.pi
        Ws = [torus_elem(3, [a, b]) for a in gr for b in gr]
    datas = [phi_data(W, lam) for W in Ws]

    mbi = np.array([mu_of(np.eye(n), dd, f, N) for dd in datas])
    ibi = int(mbi.argmin())
    print(f"  SU({N}) beta={BETA}:  bi-invariant  inf_W mu = {mbi.min():+.6f}  (at W#{ibi}, "
          f"max_W mu = {mbi.max():+.6f})")
    print(f"      formula check: -(6 c_N/N) beta + N/8 = "
          f"{-(6*(0.761882 if N==2 else 0.777245)/N)*BETA + N/8:+.6f}   [c_N from w7_pl_hess]")

    # shape family: diagonal, det 1, log-eigenvalues bounded by +-LOGB (bounded because for a
    # degenerate G the (0,2)-tensor M_G itself blows up in norm: mu -> +-inf at the boundary)
    def mu_diag(x, dd):
        return mu_of(diag_from_logx(x), dd, f, N)

    if N == 2:
        gsh = np.linspace(-2.0, 2.0, 33)
        shapes = [diag_from_logx(np.array([a, b])) for a in gsh for b in gsh]
    else:
        shapes = [np.eye(n)] + [diag_from_logx(rng.normal(0, 1.0, n - 1)) for _ in range(600)]
        for a in range(n):
            for s in (-2.0, -1.0, 1.0, 2.0):
                x = np.zeros(n - 1)
                if a < n - 1:
                    x[a] = s
                shapes.append(diag_from_logx(x))

    # per-W best over shapes
    best_per_W = np.array([min(mu_of(G, dd, f, N) for G in shapes) for dd in datas])
    # per-W max over shapes
    sup_per_W = np.array([max(mu_of(G, dd, f, N) for G in shapes) for dd in datas])
    print(f"      bounded shape family ({len(shapes)} shapes, log-eigs in [-2,2], det=1):")
    print(f"        sup_G mu  : min over W = {sup_per_W.min():+.6f}   <-- DECISIVE if < 0")
    print(f"        inf_G mu  : min over W = {best_per_W.min():+.6f}")
    print(f"        bi-invariant mu min   = {mbi.min():+.6f}")
    j = int(sup_per_W.argmin())
    print(f"        worst W for every shape: #{j};  mu(bi-inv at that W) = {mbi[j]:+.6f}")
    # boundary behaviour at the worst W
    Ww = datas[j]
    print(f"        boundary check at that W (condition number scan, SU({N})):")
    for s in (0.0, 0.5, 1.0, 2.0, 3.0, 5.0):
        x = np.zeros(n - 1); x[0] = s
        x[-1] = -s
        print(f"          log-spread {s:4.1f}: mu = {mu_of(diag_from_logx(x), Ww, f, N):+12.6f}")
    res[f"q2_SU{N}"] = dict(biinv_min=float(mbi.min()), supG_min=float(sup_per_W.min()),
                            infG_min=float(best_per_W.min()), nshapes=len(shapes))

print()
print("=" * 104)
print("STEP 5  the explicit worst configuration (realizable on the lattice)")
print("  U_1(x)=(-1)^(x2+x3+x4),  U_mu(x)=(-1)^x1 (mu=2,3,4)  on an even torus.")
print("  Then every (1,mu) plaque = -1 and every (mu,nu) plaque = +1, so a direction-1 link")
print("  sits in 6 plaques all equal to the centre element -1: Hess S_W = -(3 beta/N) delta.")
print("  The negative block is ISOTROPIC, so no metric shape can align anything with it;")
print("  (star) at t=0 then needs  Ric_G > (6 beta/N) delta.")
for N in (2, 3):
    c = 3 * BETA / N
    need, have = 2 * c, N / 4
    print(f"  SU({N}) beta={BETA}:  block = -{c:.4f} delta ;  needed lam_min(Ric_G) > {need:.4f}"
          f" ;  sup_G = {have:.4f}  ==> short by {need/have:.2f}x  (Ric must be {need/have:.2f}x larger)")
    res[f"step5_SU{N}"] = dict(block=c, need=need, have=have, factor=need / have)

with open("mc_ricci.json", "w") as fh:
    json.dump(res, fh, indent=1)
print("\nwrote mc_ricci.json")
