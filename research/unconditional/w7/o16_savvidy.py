"""O16 — is the Savvidy negative mode a friend or an enemy of H1's Var?

Object: single SU(2) plaquette, convention = research/v3/species/w7/na_plaquette.py and
research/unconditional/w7/h1e_var.py:
    U = exp(i xi.sig/2), |xi|<2pi, h=rho/2, u=xi3/rho,
    S = beta(1 - cos h cos(th/2) + sin h sin(th/2) u),  J = 4 sin^2 h / rho^2,  Phi = -log J,
    v = e3, conditionalise <v,xi> = a0 = 0  =>  xi3 = 0 (the slice).
At th=pi: S|slice == beta, G := <e3,grad S> = beta sin(rho/2)/rho, Var(G) = beta^2 V0.

Task: decompose Var(G + Phi') in the eigenbasis of Hess S_W into negative- and
positive-mode parts, get the shares at beta=2.5 and beta=6, and decide (a)/(b)/(c).

Everything below is deterministic Gauss-Legendre quadrature + closed forms.  No sampling,
no fitting, no seeds.
"""
import numpy as np, math

PI = math.pi
EULER = 0.5772156649015328606
E3 = np.array([0.0, 0.0, 1.0])


def gl(f, a, b, n=4000):
    x, w = np.polynomial.legendre.leggauss(n)
    t = 0.5 * (b - a) * x + 0.5 * (a + b)
    return 0.5 * (b - a) * float(np.sum(w * f(t)))


def Ci_num(x):
    return EULER + math.log(x) + gl(lambda t: (np.cos(t) - 1.0) / t, 0.0, x, 600)


def Cin_num(x):
    return EULER + math.log(x) - Ci_num(x)


# ------------------------------------------------------------------ [1] 1D slice constants
Z = gl(lambda x: np.sin(x) ** 2 / x, 0.0, PI)
N1 = gl(lambda x: np.sin(x) ** 3 / x ** 2, 0.0, PI)
N2 = gl(lambda x: np.sin(x) ** 4 / x ** 3, 0.0, PI)
Zc = 0.5 * Cin_num(2 * PI)
N1c = 0.75 * (math.log(3.0) + Ci_num(PI) - Ci_num(3 * PI))
N2c = math.log(2.0) + Ci_num(2 * PI) - Ci_num(4 * PI)
Ef = N1 / (2 * Z)                      # E[G]/beta
Ef2 = N2 / (4 * Z)                     # E[G^2]/beta^2
V0 = Ef2 - Ef ** 2                     # Var(G)/beta^2
V0c = (N2c * Zc - N1c ** 2) / (4 * Zc ** 2)
W0 = 1.0 / (2.0 * Cin_num(2 * PI))
BSTAR = math.sqrt(W0 / V0)

# mean curvature magnitude on the slice, |lambda_-| = beta*mu(rho), rho=2x
#   H_{i3} = (beta/rho)[cos h/2 - sin h/rho] xi_i/rho  =>  |M| = beta|x cos x - sin x|/(4x^2)
MU = gl(lambda x: np.sin(x) - x * np.cos(x), 0.0, PI, 6000) / (4.0 * Z)   # <mu>_Haar
MU = gl(lambda x: (np.sin(x) - x * np.cos(x)) * np.sin(x) ** 2 / x / (4 * x ** 2), 0.0, PI, 6000) / Z

print("=" * 96)
print("[1] slice integrals  (x = rho/2 in (0,pi); weight sin^2 x / x)")
print(f"    Z  = {Z:.15f}   closed 0.5*Cin(2pi)  = {Zc:.15f}   d = {Z-Zc:+.2e}")
print(f"    N1 = {N1:.15f}   closed = {N1c:.15f}   d = {N1-N1c:+.2e}")
print(f"    N2 = {N2:.15f}   closed = {N2c:.15f}   d = {N2-N2c:+.2e}")
print(f"    E[G]/beta      = {Ef:.15f}")
print(f"    E[G^2]/beta^2  = {Ef2:.15f}")
print(f"    V0 = Var/beta^2 = {V0:.15f}   (closed {V0c:.15f}, d={V0-V0c:+.2e})")
print(f"    W0 = 1/(2Cin 2pi) = {W0:.12f}   beta* = sqrt(W0/V0) = {BSTAR:.8f}   (H1E: 4.30660785)")
print(f"    <mu>_Haar = <|lambda_-|>/beta = {MU:.12f}")

# ------------------------------------------------------------------ geometry of S
def grad_hess_S(xi, beta, theta):
    """xi: (M,3).  Returns S (M,), grad (M,3), Hess (M,3,3).  Exact chain rule in (rho,u)."""
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
    H = H + s_u[:, None, None] * d2u
    return S, grad, H


def d3_Phi(xi):
    """d_3 Phi on/off the slice.  Phi=-log J, J=4 sin^2 h/rho^2 => d_3 Phi = -(cot h - 1/h) u."""
    xi = np.atleast_2d(np.asarray(xi, float))
    r = np.linalg.norm(xi, axis=1)
    h = 0.5 * r
    u = xi[:, 2] / r
    return -(1.0 / np.tan(h) - 1.0 / h) * u


# ------------------------------------------------------------------ [2] Hess eigen-structure on the slice
print("\n[2] Hess S_W on the conditionalisation slice {xi3=0} at theta=pi, v=e3")
rs = np.linspace(0.05, 2 * PI - 0.05, 24)
pts = np.stack([rs * np.cos(0.37 + 0.9 * np.arange(len(rs))), rs * np.sin(0.37 + 0.9 * np.arange(len(rs))),
                np.zeros_like(rs)], axis=1)
S0, G0, H0 = grad_hess_S(pts, 1.0, PI)
ev = np.linalg.eigvalsh(H0)
Mu_theory = np.abs(rs * np.cos(rs / 2) / 2.0 - np.sin(rs / 2)) / rs ** 2
print(f"    {'rho':>8}{'lam_-':>13}{'lam_0':>13}{'lam_+':>13}{'-|M| theory':>15}{'H33':>12}")
for i in range(0, len(rs), 4):
    print(f"    {rs[i]:8.4f}{ev[i,0]:13.6f}{ev[i,1]:13.6f}{ev[i,2]:13.6f}"
          f"{-Mu_theory[i]:15.6f}{H0[i,2,2]:12.1e}")
print(f"    max|H_{'{'}33{'}'}| on slice            = {np.max(np.abs(H0[:,2,2])):.3e}   (exact 0)")
print(f"    max|H_11|,|H_22|,|H_12| on slice  = {np.max(np.abs(H0[:,0,0])):.3e}, "
      f"{np.max(np.abs(H0[:,1,1])):.3e}, {np.max(np.abs(H0[:,0,1])):.3e}")
print(f"    max|lam_- + |M|| (sign+value)      = {np.max(np.abs(ev[:,0] + Mu_theory)):.3e}")
print(f"    max|lam_0|                          = {np.max(np.abs(ev[:,1])):.3e}   (azimuthal, exact 0)")
print(f"    max|lam_+ - |M||                    = {np.max(np.abs(ev[:,2] - Mu_theory)):.3e}")

# eigenvector alignment with e3
w, V = np.linalg.eigh(H0)
c2 = V[:, 2, :] ** 2                     # (M,3);  c_k = <e3, e_k> = V[m,2,k]  (eigh: column k = e_k)
print(f"    c^2 for (lam_-, lam_0, lam_+): mean = {c2.mean(axis=0)[0]:.12f}, "
      f"{c2.mean(axis=0)[1]:.3e}, {c2.mean(axis=0)[2]:.12f}   (theory 0.5, 0, 0.5)")
print(f"    max|c^2 - (1/2,0,1/2)|            = {np.max(np.abs(c2 - np.array([0.5,0.0,0.5]))):.3e}")

# ------------------------------------------------------------------ [3] exact Var + signed mode budget
def slice_moments(beta, theta, ns=700, nph=96):
    """2D quadrature on the slice {xi3=0}:  d nu ~ e^{-S} J s ds dphi.  Returns Var, mode data."""
    R = 2 * PI
    xs, ws = np.polynomial.legendre.leggauss(ns)
    xp, wp = np.polynomial.legendre.leggauss(nph)
    s = 0.5 * R * (xs + 1.0)
    wsv = 0.5 * R * ws
    ph = PI * (xp + 1.0)
    wpv = PI * wp
    S_, P_ = np.meshgrid(s, ph, indexing="ij")
    WS, WP = np.meshgrid(wsv, wpv, indexing="ij")
    xi = np.stack([(S_ * np.cos(P_)).ravel(), (S_ * np.sin(P_)).ravel(), np.zeros(S_.size)], axis=1)
    wt = (WS * S_ * WP).ravel()
    Sv, grad, H = grad_hess_S(xi, beta, theta)
    r = np.linalg.norm(xi, axis=1)
    h = 0.5 * r
    J = 4.0 * np.sin(h) ** 2 / r ** 2
    lw = np.log(wt) - Sv + np.log(J)
    mx = lw.max()
    ww = np.exp(lw - mx)
    Zs = ww.sum()
    E = lambda f: float((ww * f).sum() / Zs)
    G = grad[:, 2] + d3_Phi(xi)
    Gs = grad[:, 2]
    EG = E(G)
    Var = E(G * G) - EG ** 2
    EGs = E(Gs)
    VarG = E(Gs * Gs) - EGs ** 2
    ev, V = np.linalg.eigh(H)
    c = V[:, 2, :]                              # c_k = <e3, e_k>  (column k = eigenvector k)
    g = np.einsum("mak,ma->mk", V, grad)        # g_k = <e_k, grad S>;  V[m,:,k] = e_k

    EGk = np.array([E(g[:, k]) for k in range(3)])
    Cov = np.zeros((3, 3))
    for k in range(3):
        for l in range(3):
            Cov[k, l] = E(g[:, k] * g[:, l]) - EGk[k] * EGk[l]
    # exact split:  A_k := c_k g_k,  Var = sum_k C_k,  C_k := sum_l Cov(A_k, A_l)
    A = c * g
    Gt = A.sum(axis=1)
    EGA = np.array([E(A[:, k]) for k in range(3)])
    C = np.array([E(A[:, k] * Gt) - EGA[k] * EG for k in range(3)])
    lam = np.array([E(ev[:, k]) for k in range(3)])
    c2m = np.array([E(c[:, k] ** 2) for k in range(3)])
    return dict(Var=Var, EG=EG, C=C, lam=lam, c2=c2m, Cov=Cov, ES=E(Sv),
                dev=float(np.abs(C.sum() - Var)),
                dev_c=float(np.max(np.abs(np.array([c2m[k] * Var for k in range(3)]) - C))),
                H33=E(H[:, 2, 2]), VarG=VarG, EG2=E(G * G), phi=float(np.max(np.abs(d3_Phi(xi)))))


print("\n[3] exact Var and its Hess-eigenbasis decomposition  (Var = sum_k C_k,  C_k = c_k^2 Var)")
print(f"    {'beta':>6}{'Var(G+Phi)':>14}{'C_-(neg)':>12}{'C_0':>10}{'C_+(pos)':>12}"
      f"{'sum-Var':>11}{'frac(neg)':>11}{'frac(pos)':>11}{'Var/beta^2':>12}")
TAB = {}
for beta in (0.5, 2.5, 6.0):
    d = slice_moments(beta, PI)
    TAB[beta] = d
    fr_neg = d["C"][0] / d["Var"]
    fr_pos = d["C"][2] / d["Var"]
    print(f"    {beta:6.2f}{d['Var']:14.9f}{d['C'][0]:12.9f}{d['C'][1]:10.2e}{d['C'][2]:12.9f}"
          f"{d['dev']:11.1e}{fr_neg:11.9f}{fr_pos:11.9f}{d['Var']/beta**2:12.9f}")
print(f"    E[S] on slice at th=pi must equal beta:  " +
      ", ".join(f"{b}:{TAB[b]['ES']-b:+.1e}" for b in TAB))
print(f"    max|Phi'| on slice (exact 0)          : " + ", ".join(f"{b}:{TAB[b]['phi']:.1e}" for b in TAB))
print(f"    |Var(G+Phi')-Var(G)| on slice (exact 0): " +
      ", ".join(f"{b}:{abs(TAB[b]['Var']-TAB[b]['VarG']):.1e}" for b in TAB))
print(f"    |sum_k C_k - Var| (split residual)     : " + ", ".join(f"{b}:{TAB[b]['dev']:.1e}" for b in TAB))
print(f"    |c_k^2 Var - C_k| (C_k=c_k^2 Var check): " + ", ".join(f"{b}:{TAB[b]['dev_c']:.1e}" for b in TAB))

print("\n[4] curvature channel:  E[H] = sum_k c_k^2 lam_k   (H1E: E[H]=0 exactly at th=pi)")
print(f"    {'beta':>6}{'<lam_->':>12}{'<lam_0>':>10}{'<lam_+>':>12}{'E[H] slit':>12}"
      f"{'<H33>':>10}{'neg budget':>13}{'pos budget':>13}")
for beta in (2.5, 6.0):
    d = TAB[beta]
    EH = float(np.dot(d["c2"], d["lam"]))
    Bneg = d["c2"][0] * d["lam"][0] - d["C"][0]
    Bpos = d["c2"][2] * d["lam"][2] - d["C"][2]
    print(f"    {beta:6.2f}{d['lam'][0]:12.6f}{d['lam'][1]:10.1e}{d['lam'][2]:12.6f}{EH:12.2e}"
          f"{d['H33']:10.1e}{Bneg:13.6f}{Bpos:13.6f}")
    print(f"           theory: <|lam|>/2 = beta*<mu>/2 = {beta*MU/2:.6f} ; beta^2 V0/2 = {beta**2*V0/2:.6f}"
          f" ;  Bneg+Bpos = {Bneg+Bpos:.6f}  (-beta^2 V0 = {-beta**2*V0:.6f})")

# ------------------------------------------------------------------ [5] beta scaling of the shares
print("\n[5] beta-scaling of the negative-mode share")
print(f"    {'beta':>6}{'Var=beta^2 V0':>15}{'Var/W0':>10}{'frac(neg)':>11}{'<|lam|>/RMS(G)':>16}"
      f"{'<|lam|>/Var':>13}{'delta=W0-Var':>15}")
for beta in (0.5, 1.0, 2.5, 4.30660785, 6.0, 10.0):
    V = beta ** 2 * V0
    print(f"    {beta:6.3f}{V:15.9f}{V/W0:10.5f}{0.5:11.9f}{beta*MU/(beta*math.sqrt(V0)):16.9f}"
          f"{beta*MU/V:13.6f}{W0-V:15.9f}")

# ------------------------------------------------------------------ [6] negative-eigenvalue census over the ball
print("\n[6] negative-eigenvalue census of Hess S_W over |xi|<2pi (th=pi).  Window for 'exactly one'")
rg = np.linspace(1e-3, 2 * PI - 1e-3, 40)
th_g = np.arccos(np.linspace(-1, 1, 21))
ph_g = np.linspace(0, 2 * PI, 25, endpoint=False)
R3, T3, P3 = np.meshgrid(rg, th_g, ph_g, indexing="ij")
XI = np.stack([(R3 * np.sin(T3) * np.cos(P3)).ravel(), (R3 * np.sin(T3) * np.sin(P3)).ravel(),
               (R3 * np.cos(T3)).ravel()], axis=1)
_, _, HB = grad_hess_S(XI, 1.0, PI)
evB = np.linalg.eigvalsh(HB)
nneg = (evB < -1e-9).sum(axis=1)
for k in range(4):
    print(f"    #negative eigenvalues = {k}: {np.mean(nneg == k)*100:6.2f}% of grid points")
rr = np.linalg.norm(XI, axis=1)
for lo, hi, lab in [(0.01, PI, "rho<pi"), (PI, 2 * PI, "rho>pi")]:
    sel = (rr > lo) & (rr < hi)
    print(f"    {lab}: nneg=1 in {np.mean(nneg[sel]==1)*100:5.1f}%, nneg=2 {np.mean(nneg[sel]==2)*100:5.1f}%, "
          f"nneg=3 {np.mean(nneg[sel]==3)*100:5.1f}%")
print("    -> the 'exactly one negative direction' statement is window-dependent; see doc §5.")

# ------------------------------------------------------------------ [7] theta!=pi: shares drift with beta
print("\n[7] off the exact window: theta=0.9pi (Var/beta^2 drifts, so do the shares)")
print(f"    {'beta':>6}{'theta':>9}{'Var':>13}{'Var/beta^2':>13}{'frac(neg)':>11}{'frac(pos)':>11}")
for th in (PI, 0.9 * PI):
    for beta in (2.5, 6.0):
        d = slice_moments(beta, th)
        print(f"    {beta:6.2f}{th:9.4f}{d['Var']:13.8f}{d['Var']/beta**2:13.9f}"
              f"{d['C'][0]/d['Var']:11.9f}{d['C'][2]/d['Var']:11.9f}")

print("\nClay 0%.  No QED.")
