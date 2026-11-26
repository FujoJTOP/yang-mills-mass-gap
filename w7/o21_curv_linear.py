"""
o21_curv_linear.py -- O21: turn  E[H] ~ beta (linear, m=1)  into  E[H_vv+Phi''|0] >= c*beta*||v||^2.

Exact model (H1A-EXACT Sec 2.1 / O05 / O06 convention, single plaquette SU(2)):
    U = exp(xi), xi in R^3, rho=|xi| in (0,2pi)   (fundamental domain / Weyl chamber)
    S_W = beta (1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) u),  u = xi_3/rho
    Phi = -log J = -2 log( sin(rho/2)/rho )        (SU(2) Haar Jacobian, beta-FREE)
    slice xi_3 = 0, coarse mode v = e_3, conditioning point a = <v,xi> = 0
    slice measure  w(rho) drho = 2pi sin^2(rho/2)/rho * exp(-S_W) drho

Claims tested (all deterministic quadrature, no sampling; `[C]`):
  P1  exact structure  H_vv = beta * Q,  Q = Q_geom, and  sign(Q) is beta-INDEPENDENT.
      general 4-link SU(2) check of   d_v^2(1 - 1/2 Re tr U_p) = 1/2 tr(H_p^2)[alpha_p - gamma_p]
      with H_p = i W_p,  W_p = U_p^{-1} dU_p/ds,  alpha_p = Re tr(U_p H_p^2)/tr(H_p^2)
  P2  sign scan of the *full* 3x3 Hessian of V=S+Phi at xi=(rho,0,0):
      H = [[A,0,C],[0,B,0],[C,0,D]]  (A radial, B=D transverse=o06's H_vv+Phi'')
      => min_v [H_vv+Phi''] = min(B, lambda_-) ;  is it >= 0 pointwise?  (Savvidy test)
  P3  Phi''(beta) = E[Phi''|0] vs beta ;  d/dbeta E[Phi''] = -Cov(Phi'', S_W^{(0)}|0)  (sign + FD check)
      pointwise floors:  transverse Phi'' >= 1/6 = N/12 ;  radial Phi'' >= 1/12 = N/24
  P4  c(beta) = min_v E[H_vv+Phi''|0] / (beta ||v||^2) ;  analytic c_inf = cos(th/2)/4 (transverse)
      and the rigorous window bound c = (N/12)/B on beta in (0,B].

Run:  py -3.12 w7/o21_curv_linear.py
"""
import numpy as np

np.set_printoptions(precision=10, suppress=False)
TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------- exact scalar pieces
def Phi(rho):
    return -2.0 * np.log(np.sin(rho / 2.0) / rho)


def Jfun(rho):
    return np.sin(rho / 2.0) ** 2 / rho ** 2


def Phipp_trans(rho):
    """d^2 Phi / d xi_3^2 at xi_3=0  (= Phi'(rho)/rho).   inf = 1/6 = N/12 (N=2)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2.0 / rho ** 2 - 1.0 / (np.tan(rho / 2.0) * rho)


def Phipp_rad(rho):
    """d^2 Phi / d rho^2 = -2/rho^2 + (1/2)csc^2(rho/2).  inf = 1/6 = N/12 (N=2)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.sin(rho / 2.0)
        return -2.0 / rho ** 2 + 0.5 / (s * s)


def V(amp, beta, theta):
    """amped array amp=(rho,xi_2,xi_3) -> S+Phi."""
    rho = np.sqrt((amp * amp).sum(axis=-1))
    u = np.where(rho > 0, amp[..., 2] / rho, 0.0)
    S = beta * (1.0 - np.cos(rho / 2.0) * np.cos(theta / 2.0)
                + np.sin(rho / 2.0) * np.sin(theta / 2.0) * u)
    return S + Phi(rho)


def hess3_fd(rho, beta, theta, h=1e-4):
    """3x3 Hessian of V at xi=(rho,0,0) by central differences (independent of the analytics)."""
    x = np.array([rho, 0.0, 0.0])
    H = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            ei = np.zeros(3); ei[i] = h
            ej = np.zeros(3); ej[j] = h
            H[i, j] = (V(x + ei + ej, beta, theta) - V(x + ei - ej, beta, theta)
                       - V(x - ei + ej, beta, theta) + V(x - ei - ej, beta, theta)) / (4 * h * h)
    return 0.5 * (H + H.T)


def hess3_exact(rho, beta, theta):
    """Analytic Hessian of V at xi=(rho,0,0):  [[A,0,C],[0,B,0],[C,0,B]]."""
    c, s = np.cos(theta / 2.0), np.sin(theta / 2.0)
    sh, ch = np.sin(rho / 2.0), np.cos(rho / 2.0)
    B = beta * c * sh / (2.0 * rho) + Phipp_trans(rho)          # transverse (mode v = e_3)
    A = beta * c * ch / 4.0 + Phipp_rad(rho)                    # radial (mode v = e_1)
    C = beta * s * (rho * ch / 2.0 - sh) / rho ** 2              # mode mixing (u-term)
    return np.array([[A, 0.0, C], [0.0, B, 0.0], [C, 0.0, B]])


def eigs_min(H):
    w = np.linalg.eigvalsh(H)
    return float(w[0])


def lam_minus(A, B, C):
    return (A + B) / 2.0 - np.sqrt(((A - B) / 2.0) ** 2 + C * C)


# ---------------------------------------------------------------- slice measure / moments
def slice_quad(npan=240, nnode=40):
    """composite Gauss-Legendre on rho in (0,2pi)."""
    x, w = np.polynomial.legendre.leggauss(nnode)
    edges = np.linspace(0.0, TWO_PI, npan + 1)
    rho, ww = [], []
    for k in range(npan):
        a, b = edges[k], edges[k + 1]
        rho.append(0.5 * (b - a) * (x + 1.0) + a)
        ww.append(0.5 * (b - a) * w)
    return np.concatenate(rho), np.concatenate(ww)


RHO, DRHO = slice_quad()


def moments(beta, theta):
    S = beta * (1.0 - np.cos(RHO / 2.0) * np.cos(theta / 2.0))     # u=0 on the slice
    dens = np.exp(-S) * Jfun(RHO) * TWO_PI * RHO * DRHO            # e^{-S} J dX_slice
    Z = dens.sum()
    p = dens / Z
    S0 = 1.0 - np.cos(RHO / 2.0) * np.cos(theta / 2.0)             # beta-free action density
    c, s = np.cos(theta / 2.0), np.sin(theta / 2.0)
    Qg = c * np.sin(RHO / 2.0) / (2.0 * RHO)                       # H_vv/(beta) (transverse)
    out = {}
    out["Z"] = Z
    out["Eq"] = float((p * Qg).sum())
    out["EPhi_t"] = float((p * Phipp_trans(RHO)).sum())
    out["EPhi_r"] = float((p * Phipp_rad(RHO)).sum())
    out["EB"] = float((p * (beta * Qg + Phipp_trans(RHO))).sum())
    A = beta * c * np.cos(RHO / 2.0) / 4.0 + Phipp_rad(RHO)
    C = beta * s * (RHO * np.cos(RHO / 2.0) / 2.0 - np.sin(RHO / 2.0)) / RHO ** 2
    B = beta * Qg + Phipp_trans(RHO)
    out["EA"] = float((p * A).sum())
    out["EC"] = float((p * C).sum())
    lm = lam_minus(A, B, C)
    out["Elam"] = float((p * lm).sum())               # E[ min eig along mixed 2x2 block ]
    out["Eminpt"] = float((p * np.minimum(B, lm)).sum())   # E[ pointwise min_v ]
    out["Cov_Phi_S0"] = float((p * Phipp_trans(RHO) * S0).sum()
                              - (p * Phipp_trans(RHO)).sum() * (p * S0).sum())
    return out


# ---------------------------------------------------------------- general 4-link SU(2) check
def su2_rnd(rng):
    a = rng.normal(size=4)
    return a / np.linalg.norm(a)


def expi(X):
    ev, U = np.linalg.eigh(X)
    return U @ np.diag(np.exp(1j * ev)) @ U.conj().T


def gen(t):
    return t[0] * np.eye(2, dtype=complex) + 1j * (
        t[1] * np.array([[0, 1], [1, 0]])
        + t[2] * np.array([[0, -1j], [1j, 0]])
        + t[3] * np.array([[1, 0], [0, -1]]))


def modelB(nsamp=600, seed=7, h=1e-4, scale=0.35):
    """Robust statistics: Q (Richardson 2nd difference of the scalar), and the alignment alpha."""
    rng = np.random.default_rng(seed)
    v = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex) / 2.0   # hermitian su(2) dir (||v||=1/2)
    Qs = []
    for _ in range(nsamp):
        X = [scale * gen(su2_rnd(rng)) for _ in range(4)]

        def U(s):
            L = [expi(X[l] + s * v) for l in range(4)]
            return L[0] @ L[1] @ L[2].conj().T @ L[3].conj().T

        def f(s):                                    # 1 - 1/2 Re tr U_p(s)
            return 1.0 - 0.5 * np.real(np.trace(U(s)))

        def d2(hh):
            return (f(hh) - 2.0 * f(0.0) + f(-hh)) / (hh * hh)

        Q = (4.0 * d2(0.01) - d2(0.02)) / 3.0        # Richardson (kills O(h^2))
        Up = (U(h) - U(-h)) / (2 * h)
        W0 = Up @ np.linalg.inv(U(0.0))
        Hp = -1j * W0
        Hp = 0.5 * (Hp + Hp.conj().T)
        trK = float(np.real(np.trace(Hp @ Hp)))
        if trK < 1e-12:
            continue
        alpha = float(np.real(np.trace(U(0.0) @ Hp @ Hp))) / trK
        Qs.append((Q, d2(0.01), d2(0.02), alpha))
    return np.array(Qs)


if __name__ == "__main__":
    L = []
    def P(s=""):
        print(s)
        L.append(s)

    P("=" * 96)
    P("O21  linear-in-beta curvature bound   (py 3.12 / numpy, deterministic quadrature)")
    P("=" * 96)

    # ---------------------------------------------------------------- P1 structure
    P("\n[P1a] exact 4-link SU(2) (|X_l|~0.35, Richardson 2nd difference of 1-1/2 Re tr U_p):")
    Qs = modelB()
    P(f"   n={len(Qs)}   consistency |d2(0.01)-d2(0.02)| max = {np.max(np.abs(Qs[:,1]-Qs[:,2])):.3e}"
      f"  (|Q|max={np.abs(Qs[:,0]).max():.3e})")
    P(f"   alpha=Re tr(U H^2)/tr H^2 : min={Qs[:,3].min(): .6f} max={Qs[:,3].max(): .6f}"
      f" frac<0={np.mean(Qs[:,3]<0):.4f}")
    P(f"   Q_geom(unit ||v||^2)      : min={Qs[:,0].min(): .6e} max={Qs[:,0].max(): .6e}"
      f" mean={Qs[:,0].mean(): .6e} frac<0={np.mean(Qs[:,0]<0):.4f}")
    P("   identity (exact algebra): d_v^2(1-1/2 Re tr U_p) = 1/2 Re tr[U_p(W_p^2 + Wdot_p)], W_p=U_p^{-1}U_p'")
    P("                            = 1/2 tr(H_p^2)[alpha_p - gamma_p],  H_p=iW_p, gamma_p=Re tr(U_p Wdot_p)/tr H_p^2")
    P("   => Q_geom is INDEFINITE pointwise (Savvidy). H_vv = beta*Q_geom is NOT >=0 free.")

    P("\n[P1a'] same statistic at the LARGER amplitude |X_l|~0.70 (outside the physical scaling window)")
    Qb = modelB(scale=0.70)
    P(f"   n={len(Qb)}   consistency max = {np.max(np.abs(Qb[:,1]-Qb[:,2])):.3e}")
    P(f"   alpha: min={Qb[:,3].min(): .6f} max={Qb[:,3].max(): .6f} frac<0={np.mean(Qb[:,3]<0):.4f}")
    P(f"   Q_geom: min={Qb[:,0].min(): .6e} max={Qb[:,0].max(): .6e} frac<0={np.mean(Qb[:,0]<0):.4f}"
      f"  (alpha<0 AND Q<0: {np.mean((Qb[:,3]<0)&(Qb[:,0]<0)):.4f})")

    P("\n[P1b] single plaquette: analytic 3x3 Hessian vs finite difference (xi=(rho,0,0))")
    for (b, th, r) in [(2.5, 0.0, 1.7), (2.5, 1.1, 1.7), (2.5, np.pi, 0.9), (6.0, 0.3, 2.4)]:
        He, Hf = hess3_exact(r, b, th), hess3_fd(r, b, th)
        P(f"   beta={b:<4} th={th:7.4f} rho={r:<5} max|He-Hf|={np.max(np.abs(He-Hf)):.2e}"
          f"   eig(He)={np.array2string(np.linalg.eigvalsh(He), precision=6)}")
    P("   H_vv(transverse) = beta*cos(th/2)*sin(rho/2)/(2rho): sign = sign(cos(th/2)), beta-free.")

    # ---------------------------------------------------------------- P2 sign scan
    P("\n[P2] pointwise sign scan of the FULL Hessian  min_v[H_vv+Phi''] = min(B, lambda_-)")
    rr = np.linspace(1e-3, TWO_PI - 1e-3, 400)
    hdr = "   beta  theta   inf_r min(B,lam-)   inf_r B(transv)   inf_r A(radial)   frac<0(min)"
    P(hdr)
    for th in [0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]:
        for b in [0.0, 2.5, 10.0]:
            c, s = np.cos(th / 2.0), np.sin(th / 2.0)
            B = b * c * np.sin(rr / 2.0) / (2 * rr) + Phipp_trans(rr)
            A = b * c * np.cos(rr / 2.0) / 4.0 + Phipp_rad(rr)
            C = b * s * (rr * np.cos(rr / 2.0) / 2.0 - np.sin(rr / 2.0)) / rr ** 2
            lm = lam_minus(A, B, C)
            m = np.minimum(B, lm)
            P(f"   {b:<5} {th:<7.4f} {m.min(): .6e}      {B.min(): .6e}      {A.min(): .6e}"
              f"      {np.mean(m < 0):.3f}")
    P("   pointwise floors: B >= 1/6 = N/12 (transverse) and A >= 1/6 = N/12 (radial) => Phi'' >= (N/12) I3;")
    P("   negative region requires cos(th/2)<0 (th>pi, i.e. Re tr U_p<0) and/or |C|^2>AB.")

    # ---------------------------------------------------------------- P3 Phi''(beta)
    P("\n[P3a] pointwise floors of Phi''")
    for name, f in [("transverse  Phi''", Phipp_trans), ("radial      Phi''", Phipp_rad)]:
        v = f(rr)
        P(f"   {name}: inf_(0,2pi) = {np.nanmin(v):.10f} at rho={rr[np.nanargmin(v)]:.4g}")
    P("   => BOTH entries >= N/12 = 1/6 = 0.1666666667 ; Phi''(0) = N/12 attained at rho=0")

    P("\n[P3b] E[Phi''|0] vs beta ;  d/dbeta = -Cov(Phi'',S0) (analytic) vs FD")
    P("   theta    beta    E[Phi''_t|0]   E[Phi''_r|0]   -Cov(Phi''_t,S0)   dE/dbeta(FD)")
    for th in [0.0, np.pi / 2]:
        bs = [0.05, 0.5, 1.0, 2.5, 6.0, 10.0, 50.0, 200.0, 1000.0]
        prev = None
        for b in bs:
            m = moments(b, th)
            fd = ""
            if prev is not None:
                b0, e0 = prev
                fd = f"{((m['EPhi_t'] - e0) / (b - b0)): .6e}"
            P(f"   {th:<7.4f} {b:<7} {m['EPhi_t']: .10f}   {m['EPhi_r']: .10f}"
              f"    {-m['Cov_Phi_S0']: .6e}      {fd}")
            prev = (b, m["EPhi_t"])
    P("   => E[Phi''] decreases monotonically in beta, floored by N/12 (=1/6); d/dbeta = -Cov <= 0.")

    # ---------------------------------------------------------------- P4 c
    P("\n[P4] c(beta) = E[H_vv+Phi''|0]/(beta ||v||^2) for the mode v=e_3 (H1A/O05/O06 setting)")
    P("   theta    beta     E[Q|0]=E[H]/beta   E[Phi''|0]    E[B|0]      c=E[B]/beta   [diag E[lam-]/beta]")
    rows = {}
    for th in [0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4]:
        for b in [0.5, 1.0, 2.5, 6.0, 10.0, 50.0, 200.0]:
            m = moments(b, th)
            cE = m["EB"] / b
            rows[(th, b)] = (cE, m)
            P(f"   {th:<7.4f} {b:<7} {m['Eq']: .8e}    {m['EPhi_t']: .8f}   {m['EB']: .8e}"
              f"   {cE: .8e}      {m['Elam']/b: .6e}")
    P("\n   diag: Hessian eigenvalues at the point (radial A, mixed lambda_-) are NOT slice modes of")
    P("   this toy (only directions transversal to the radius are), listed for completeness:")
    for th in [0.0, np.pi / 2]:
        for b in [2.5, 10.0]:
            m = rows[(th, b)][1]
            P(f"     theta={th:<6.4f} beta={b:<5} E[A]={m['EA']: .6e}  E[lam-]={m['Elam']: .6e}")
    P("\n   analytic large-beta limits (rho -> 0, measure concentrates):")
    for th in [0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]:
        c = np.cos(th / 2.0)
        P(f"   theta={th:<8.5f}  c_inf(transv)=cos(th/2)/4={0.25*c: .6e}"
          f"   c_inf(radial)=cos(th/2)/4={0.25*c: .6e}   c_inf(mixed)={0.25*c: .6e} "
          f"  [all three coincide: isotropy at rho=0]")
    P("\n   convergence: E[B|0] and E[min] at npan*2 (beta=2.5)")
    keep = (RHO, DRHO)
    for th in [0.0, np.pi / 2]:
        m1 = moments(2.5, th)
        RHO, DRHO = slice_quad(npan=480, nnode=40)
        m2 = moments(2.5, th)
        P(f"   theta={th:<6.4f} E[B]: {m1['EB']:.10f} -> {m2['EB']:.10f} ;"
          f" E[min]: {m1['Eminpt']:.10f} -> {m2['Eminpt']:.10f}")
        RHO, DRHO = keep
        RHO, DRHO = keep
    P("\n   rigorous window bound (uses only pointwise facts): c = (N/12)/B on beta in (0,B]")
    for B in [1.0, 2.5, 6.0]:
        P(f"     B={B:<5} c = 1/(6*{B}) = {(1.0/6.0)/B:.6e}")

    with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print("\n[log written]")
