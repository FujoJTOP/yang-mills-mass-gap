"""Q12-SOS -- can  S - D = E[H+Phi''] - Var(G+Phi')  be written as a sum of nonnegative terms?

Conventions copied from research/unconditional/w7/o16_savvidy.py (single SU(2) plaquette):
    U = exp(i xi.sig/2), S = beta(1 - cos h cos(th/2) + sin h sin(th/2) u), h=rho/2, u=xi3/rho,
    J = 4 sin^2 h / rho^2, Phi = -log J, v = e3 (unit: <v,v> = 2tr T3^2 = 1), N/12 = 1/6.
Conditioning: <v,xi> = a  =>  xi3 = a   (a=0 is o16's slice; a!=0 is the P21/P22 slice).

IMPORTANT (found here): o16's *analytic* grad/Hess of S is exact on the slice xi3=0 (all o16 checks
pass there) but WRONG off it -- FD residual rel 1.9-3.3 at xi3!=0, see [0].  So this script computes
grad S and Hess S by central differences of S (h=1e-5 / 1e-4, truncation O(h^2)<=1e-8).  Phi is
still closed-form (Phi_rho = 2/rho - cot(rho/2), Phi_rr = -2/rho^2 + csc^2(rho/2)/2).
Deterministic Gauss-Legendre quadrature, no sampling, no seeds, no fitting.
"""
import numpy as np, math

PI = math.pi
N12 = 1.0 / 6.0
E3 = np.array([0.0, 0.0, 1.0])


def S_only(xi, beta, theta):
    xi = np.atleast_2d(np.asarray(xi, float))
    r = np.linalg.norm(xi, axis=1)
    u = xi[:, 2] / r
    h = 0.5 * r
    c = math.cos(theta / 2.0)
    su = math.sin(theta / 2.0)
    return beta * (1.0 - np.cos(h) * c + np.sin(h) * su * u)


def grad_hess_S_ana(xi, beta, theta):
    """o16's analytic formulas, kept verbatim for the [0] check only."""
    xi = np.atleast_2d(np.asarray(xi, float))
    M = xi.shape[0]
    r = np.linalg.norm(xi, axis=1)
    u = xi[:, 2] / r
    h = 0.5 * r
    c = math.cos(theta / 2.0)
    su = math.sin(theta / 2.0)
    S = beta * (1.0 - np.cos(h) * c + np.sin(h) * su * u)
    s_r = beta * (0.5 * np.sin(h) * c + 0.5 * np.cos(h) * su * u)
    s_u = beta * np.sin(h) * su
    s_rr = beta * (0.25 * np.cos(h) * c - 0.25 * np.sin(h) * su * u)
    s_ru = beta * 0.5 * np.cos(h) * su
    rhat = xi / r[:, None]
    gu = (np.broadcast_to(E3, xi.shape) - u[:, None] * rhat) / r[:, None]
    grad = s_r[:, None] * rhat + s_u[:, None] * gu
    H = (s_rr[:, None, None] * rhat[:, :, None] * rhat[:, None, :]
         + (s_r / r)[:, None, None] * (np.eye(3)[None] - rhat[:, :, None] * rhat[:, None, :])
         + s_ru[:, None, None] * (rhat[:, :, None] * gu[:, None, :] + gu[:, :, None] * rhat[:, None, :]))
    xi3 = xi[:, 2]
    d2u = np.zeros((M, 3, 3))
    for a in range(3):
        for b in range(3):
            d2u[:, a, b] = (-((1.0 if a == 2 else 0.0) * xi[:, b] + xi3 * (1.0 if a == b else 0.0)) / r ** 3
                            + 3.0 * xi3 * xi[:, a] * xi[:, b] / r ** 5)
    return S, grad, H + s_u[:, None, None] * d2u


def grad_S_fd(xi, beta, theta, h=1e-5):
    g = np.zeros((xi.shape[0], 3))
    for k in range(3):
        d = np.zeros(3)
        d[k] = h
        g[:, k] = (S_only(xi + d, beta, theta) - S_only(xi - d, beta, theta)) / (2 * h)
    return g


def hess_S_fd(xi, beta, theta, h=1e-4):
    n = xi.shape[0]
    H = np.zeros((n, 3, 3))
    for a in range(3):
        for b in range(a, 3):
            da, db = np.zeros(3), np.zeros(3)
            da[a] = h
            db[b] = h
            if a == b:
                v = (S_only(xi + da, beta, theta) - 2 * S_only(xi, beta, theta)
                     + S_only(xi - da, beta, theta)) / h ** 2
            else:
                v = (S_only(xi + da + db, beta, theta) - S_only(xi + da - db, beta, theta)
                     - S_only(xi - da + db, beta, theta) + S_only(xi - da - db, beta, theta)) / (4 * h ** 2)
            H[:, a, b] = v
            H[:, b, a] = v
    return H


def Phi_d3(xi):
    xi = np.atleast_2d(np.asarray(xi, float))
    r = np.linalg.norm(xi, axis=1)
    return (2.0 / r - 1.0 / np.tan(0.5 * r)) * xi[:, 2] / r


def Phi_vv(xi):
    xi = np.atleast_2d(np.asarray(xi, float))
    r = np.linalg.norm(xi, axis=1)
    u = xi[:, 2] / r
    fp = 2.0 / r - 1.0 / np.tan(0.5 * r)
    fpp = -2.0 / r ** 2 + 0.5 / np.sin(0.5 * r) ** 2
    return fpp * u ** 2 + (fp / r) * (1.0 - u ** 2)


def jac_Phi(xi):
    xi = np.atleast_2d(np.asarray(xi, float))
    return 4.0 * np.sin(0.5 * np.linalg.norm(xi, axis=1)) ** 2 / np.linalg.norm(xi, axis=1) ** 2


def slice_pts(a, ns=700, nph=96):
    smax = math.sqrt(max((2 * PI) ** 2 - a * a, 1e-12))
    xs, ws = np.polynomial.legendre.leggauss(ns)
    xp, wp = np.polynomial.legendre.leggauss(nph)
    s = 0.5 * smax * (xs + 1.0)
    wsv = 0.5 * smax * ws
    ph = PI * (xp + 1.0)
    SS, PP = np.meshgrid(s, ph, indexing="ij")
    WS, WP = np.meshgrid(wsv, PI * wp, indexing="ij")
    xi = np.stack([(SS * np.cos(PP)).ravel(), (SS * np.sin(PP)).ravel(), np.full(SS.size, a)], axis=1)
    return xi, (WS * SS * WP).ravel()


def slice_moments(a, beta, theta, fd_hess=True):
    xi, wt = slice_pts(a)
    Sv = S_only(xi, beta, theta)
    grad = grad_S_fd(xi, beta, theta)
    H = hess_S_fd(xi, beta, theta) if fd_hess else grad_hess_S_ana(xi, beta, theta)[2]
    J = jac_Phi(xi)
    lw = np.log(wt) - Sv + np.log(J)
    ww = np.exp(lw - lw.max())
    Z = ww.sum()
    E = lambda f: float((ww * f).sum() / Z)

    G = grad[:, 2] + Phi_d3(xi)
    Gs = grad[:, 2]
    Hvv = H[:, 2, 2]
    Ph2 = Phi_vv(xi)
    EG, EGs = E(G), E(Gs)
    EPh1 = E(Phi_d3(xi))
    VarF = E(G * G) - EG ** 2
    VarG = E(Gs * Gs) - EGs ** 2
    VarPh1 = E(Phi_d3(xi) ** 2) - EPh1 ** 2
    Cov = E(Gs * Phi_d3(xi)) - EGs * EPh1

    ev, V = np.linalg.eigh(H)
    c = V[:, 2, :]
    g = np.einsum("mak,ma->mk", V, grad)
    A = c * g
    Gt = A.sum(axis=1)
    EGA = np.array([E(A[:, k]) for k in range(3)])
    C = np.array([E(A[:, k] * Gt) - EGA[k] * EG for k in range(3)])
    lam = np.array([E(ev[:, k]) for k in range(3)])
    c2 = np.array([E(c[:, k] ** 2) for k in range(3)])
    return dict(a=a, VarF=VarF, VarG=VarG, VarPh1=VarPh1, Cov=Cov, EG=EG,
                EHvv=E(Hvv), EPh2=E(Ph2), C=C, lam=lam, c2=c2,
                dev=float(abs(C.sum() - VarF)), devc=float(np.max(np.abs(c2 * VarF - C))),
                fracNeg=float(np.mean(ev[:, 0] < 0)),
                lam_pt_neg=float(np.min(ev[:, 0])))


def mG_deriv(a, beta, theta, h=2e-4):
    return (slice_moments(a + h, beta, theta)["EG"]
            - slice_moments(a - h, beta, theta)["EG"]) / (2 * h)


print("=" * 100)
print("Q12-SOS | S-D = E[H+Phi''] - Var(G+Phi') | single SU(2) plaquette, v=e3, <v,v>=1, N/12=1/6")
print("=" * 100)

# ------------------------------------------------------------------ [0] is o16's analytic Hess OK?
print("\n[0] o16 analytic grad/Hess vs central differences of S  (h=1e-5 / 1e-4)")
pts = np.array([[0.9, 0.4, 0.0], [0.9, 0.4, 0.75], [1.3, -0.6, 0.2], [0.1, 0.05, -1.1]])
Sa, ga, Ha = grad_hess_S_ana(pts, 2.5, PI)
gf, Hf = grad_S_fd(pts, 2.5, PI), hess_S_fd(pts, 2.5, PI)
print(f"    {'xi':>26}{'rel|grad2 ana-FD|':>20}{'rel|H33 ana-FD|':>18}{'rel|H_11 ana-FD|':>18}")
for i in range(len(pts)):
    rg = abs(ga[i, 2] - gf[i, 2]) / max(abs(gf[i, 2]), 1e-12)
    r3 = abs(Ha[i, 2, 2] - Hf[i, 2, 2]) / max(abs(Hf[i, 2, 2]), 1e-12)
    r1 = abs(Ha[i, 0, 0] - Hf[i, 0, 0]) / max(abs(Hf[i, 0, 0]), 1e-12)
    print(f"    {str(np.round(pts[i],3)):>26}{rg:20.3e}{r3:18.3e}{r1:18.3e}")
print("    -> analytic Hess of S is exact only on the xi3=0 slice; FD is used below.")

# ------------------------------------------------------------------ [1] 4-term identity check
print("\n[1] algebraic rewrite  S-D = E[H+Phi''] - E[(G+Phi')^2] + E[G+Phi']^2   (Var def, [T] algebra)")
print("    cross-check vs p21 (Q-id):  S-D = m_G' + E[Phi''] - Var(Phi') - Cov(G,Phi'),  m_G'=d/da E[G|a]")
print(f"    {'a':>6}{'beta':>6}{'S-D':>15}{'(Q-id)':>15}{'diff':>11}{'VarF':>12}{'EHvv':>12}"
      f"{'EPhi2':>11}{'mG':>12}")
ROW = {}
for a in (0.0, 0.25, 0.5, 0.75, 1.0):
    d = slice_moments(a, 2.5, PI)
    mG = mG_deriv(a, 2.5, PI)
    sd = d["EHvv"] + d["EPh2"] - d["VarF"]
    qid = mG + d["EPh2"] - d["VarPh1"] - d["Cov"]
    ROW[a] = dict(d=d, mG=mG, sd=sd, qid=qid)
    print(f"    {a:6.2f}{2.5:6.2f}{sd:15.9f}{qid:15.9f}{sd - qid:11.1e}{d['VarF']:12.9f}"
          f"{d['EHvv']:12.6f}{d['EPh2']:11.6f}{mG:12.6f}")
print("    (both sides use the same FD quadrature; diff is the residual of the bookkeeping)")
print("    P21 numbers for comparison (beta=2.5,th=pi): m_G'(0.25)=-0.13125, m_G'(0.75)=-0.23986,")
print("      Cov(0.75)=-1.1406e-2, Var(Phi')(0.75)=3.610e-3, E[Phi''](0.75)=0.226178.")

# ------------------------------------------------------------------ [2] per-mode ledger
print("\n[2] per-mode ledger in the Hess eigenbasis.  Canonical (floor-proportional) split:")
print("    Hess Phi >= (N/12) Id  =>  Phi''_k := c_k^2 N/12  and")
print("    S-D = R + sum_k c_k^2 (lam_k + N/12 - VarF),   R := E[Phi''] - N/12 >= 0 (P16).")
for (a, beta) in ((0.0, 2.5), (0.0, 6.0), (0.75, 2.5)):
    d = slice_moments(a, beta, PI)
    R = d["EPh2"] - N12
    tot = sum(d["c2"][k] * (d["lam"][k] + N12 - d["VarF"]) for k in range(3))
    print(f"    a={a:5.2f} beta={beta:4.2f}: VarF={d['VarF']:.9f} R={R:+.9f} "
          f"sum_k[...]={tot:+.9f}  sum+R={tot + R:+.9f}  (S-D={d['EHvv'] + d['EPh2'] - d['VarF']:+.9f})")
    print(f"        {'k':>3}{'c^2':>11}{'<lam_k>':>12}{'C_k':>12}{'c^2 VarF':>12}"
          f"{'c^2(lam+N/12-VarF)':>20}{'lam+N/12-VarF':>16}")
    for k in range(3):
        print(f"        {k:3d}{d['c2'][k]:11.6f}{d['lam'][k]:12.6f}{d['C'][k]:12.9f}"
              f"{d['c2'][k]*d['VarF']:12.9f}{d['c2'][k]*(d['lam'][k]+N12-d['VarF']):20.9f}"
              f"{d['lam'][k] + N12 - d['VarF']:16.9f}")
    print(f"        checks: |sum_k C_k-VarF|={d['dev']:.1e}  max|c^2 VarF-C_k|={d['devc']:.1e}  "
          f"frac(lam_-,k=0 < 0)={d['fracNeg']:.3f}  pointwise min lam_-={d['lam_pt_neg']:.4f}")

# ------------------------------------------------------------------ [3] pairing tests
print("\n[3] pairing tests for the negative mode (k=0, lam_- = -|lam_-| < 0):")
print("    T1 canonical floor split : VarF <= N/12 - <|lam_->|   (necessary AND sufficient mode-wise)")
print("    T2 free allocation       : nonneg per mode with free transfers <=> global S-D >= 0 (vacuous)")
print(f"    {'a':>6}{'beta':>6}{'VarF':>12}{'<|lam_->|>':>13}{'N/12-|lam_-|':>15}{'T1':>7}"
      f"{'S-D':>13}{'T2':>7}{'lam_k+N/12-VarF k=0':>21}")
for a in (0.0, 0.25, 0.5, 0.75, 1.0):
    d = ROW[a]["d"]
    Ab = max(-d["lam"][0], 0.0)
    print(f"    {a:6.2f}{2.5:6.2f}{d['VarF']:12.9f}{Ab:13.6f}{N12 - Ab:15.6f}"
          f"{str(d['VarF'] <= N12 - Ab):>7}{ROW[a]['sd']:13.6f}{str(ROW[a]['sd'] >= 0):>7}"
          f"{d['lam'][0] + N12 - d['VarF']:21.9f}")
print("\n    beta window of T1 on the a=0 slice (exact scaling Var=beta^2 V0 = 0.011059294 beta^2,")
print("    <|lam_->| = beta*0.081337146669):  solve  beta*<mu> + beta^2*V0 <= N/12")
mu, V0 = 0.081337146669, 0.011059293772974
disc = mu ** 2 + 4 * V0 * N12
bc = (-mu + math.sqrt(disc)) / (2 * V0)
print(f"      beta_c = (-mu+sqrt(mu^2+4 V0 N/12))/(2 V0) = {bc:.6f}")
for beta in (1.0, 1.5, bc, 2.5, 4.30660785, 6.0):
    print(f"      beta={beta:8.4f}: beta*<mu>+beta^2*V0 = {beta * mu + beta ** 2 * V0:.6f}"
          f"   vs N/12 = {N12:.6f}   -> T1 {beta * mu + beta ** 2 * V0 <= N12}")

# ------------------------------------------------------------------ [4] o06 pointwise witness
print("\n[4] o06 pointwise witness (model B, beta=2.5, sum_l ||v_l||^2 = 2): min H_vv = -1.893109")
print(f"    => min lam_- per unit ||v||^2 = {-1.893109 / 2:.6f};  ratio to the floor:")
print(f"    |min lam_-| / (N/12) = {1.893109 / 2 / N12:.4f}  > 1  => no Jacobian floor can pair it.")
print("\nClay 0%.  No QED.")
