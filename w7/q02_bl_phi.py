"""q02_bl_phi.py -- Q02: Hess Phi >= (N/12) I  -->  Brascamp-Lieb variance bound  -->  is it enough?

Single plaquette SU(2), beta=2.5, slice xi_3=0, v=e_3.  py -3.12, numpy only, deterministic quadrature.

Exact pieces (o21 convention, reused verbatim; already [T] in w7/o21_curv_linear.py):
    S      = beta(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) u),  u = xi_3/rho
    Phi    = -2 log(sin(rho/2)/rho)                      (Haar log-Jacobian, flat chart)
    G_v    = d S/d xi_3 |0 = beta sin(th/2) sin(rho/2)/rho
    H_vv   = d^2 S/d xi_3^2 |0 = beta cos(th/2) sin(rho/2)/(2 rho)
    Phi''  = d^2 Phi/d xi_3^2 |0 = 2/rho^2 - cot(rho/2)/rho        (inf = 1/6 = N/12, N=2)
    Phi'   = 0 at xi_3 = 0 (parity)

Brascamp-Lieb (J. Funct. Anal. 22 (1976) 366): mu = e^{-V}dX, Hess V >= lam I  =>
    Var_mu(f) <= E_mu<grad f, (Hess V)^{-1} grad f> <= (1/lam) E_mu|grad f|^2.
Here grad is the FULL link-space gradient; on the slice reading grad f = d_s f:
    f = G_v  =>  grad f = Hess S v,  |grad G_v|^2 = <Hess S v, Hess S v> >= H_vv^2  (Cauchy-Schwarz)
    f = Phi' =>  grad f = Hess Phi v, |grad Phi'|^2 = <Hess Phi v, Hess Phi v>

Output is written to q02_bl_phi.log by the runner.
"""
import numpy as np

TWO_PI = 2.0 * np.pi
N = 2.0
LAM = N / 12.0            # 1/6 = 0.1666667  [T] p16 [3]
INV_LAM = 1.0 / LAM       # 6 = 12/N


def Phi(rho):
    return -2.0 * np.log(np.sin(rho / 2.0) / rho)


def Jfun(rho):
    return np.sin(rho / 2.0) ** 2 / rho ** 2


def Phipp_trans(rho):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 2.0 / rho ** 2 - 1.0 / (np.tan(rho / 2.0) * rho)


def Phipp_rad(rho):
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.sin(rho / 2.0)
        return -2.0 / rho ** 2 + 0.5 / (s * s)


def slice_quad(npan=240, nnode=40):
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
    c, s = np.cos(theta / 2.0), np.sin(theta / 2.0)
    S = beta * (1.0 - np.cos(RHO / 2.0) * c)
    dens = np.exp(-S) * Jfun(RHO) * TWO_PI * RHO * DRHO
    p = dens / dens.sum()
    H = beta * c * np.sin(RHO / 2.0) / (2.0 * RHO)
    G = beta * s * np.sin(RHO / 2.0) / RHO
    P2 = Phipp_trans(RHO)
    def m(x):
        return float((p * x).sum())
    out = dict(EH=m(H), EH2=m(H * H), VarH=m((H - m(H)) ** 2),
               EG=m(G), VarG=m((G - m(G)) ** 2),
               EP2=m(P2), EP22=m(P2 * P2), VarP2=m((P2 - m(P2)) ** 2),
               EP1=m(np.zeros_like(RHO)))
    return out


OUT = []
def P(t=""):
    print(t)
    OUT.append(t)


P("=" * 96)
P("Q02-BL-PHI  Brascamp-Lieb with Hess(Phi) >= (N/12) I  ->  variance bound?   SU(2) N=2, 12/N=%.1f" % INV_LAM)
P("=" * 96)

# ---------------------------------------------------------------- [0] premise: the floor itself
P("\n[0] premise (i): Hess Phi >= (N/12) I  [T, p16 [3]]; re-check the two eigenvalues on (0,2pi)")
grid = np.linspace(1e-4, TWO_PI - 1e-4, 400001)
P("    inf Phi''_trans = %.10f   inf Phi''_rad = %.10f   N/12 = %.10f"
  % (Phipp_trans(grid).min(), Phipp_rad(grid).min(), LAM))
P("    => floor holds in BOTH directions, ratio to N/12 = %.6f / %.6f"
  % (Phipp_trans(grid).min() / LAM, Phipp_rad(grid).min() / LAM))

# ---------------------------------------------------------------- [1] premise: the interacting Hessian
P("\n[1] premise (ii): BL needs Hess V >= lam I for V = S_W + Phi (the density e^{-V} = e^{-S_W} dHaar)")
P("    Hess V = Hess S_W + Hess Phi;  Hess Phi >= 1/6 but Hess S_W is Savvidy-indefinite.")
P("    measured pointwise min H_vv = -1.893109 (p16 log [7], o06 [4], beta=2.5)")
P("    full Hessian min at beta=2.5 (lattice-cpb REPORT sec.3.4/3.5): kappa_eff < 0, min ~ -2.8")
P("    => V is NOT lambda-convex; BL hypothesis FAILS for the Wilson measure (verdict (c)).")

# ---------------------------------------------------------------- [2] slice moments + BL numbers
P("\n[2] single plaquette beta=2.5, slice xi_3=0, v=e_3: moments and BL right-hand sides")
P("    th      E[H_vv]     E[H_vv^2]   (12/N)E[H^2] | E[Phi'']   E[Phi''^2]  (12/N)E[Phi''^2]")
rows = {}
for th in (0.0, np.pi / 2, np.pi):
    m = moments(2.5, th)
    rows[th] = m
    P("   %.3f  %10.6f  %10.6f  %11.4f  | %8.5f  %10.5f  %12.4f"
      % (th, m["EH"], m["EH2"], INV_LAM * m["EH2"], m["EP2"], m["EP22"], INV_LAM * m["EP22"]))
P("    (slice reading: |grad G_v|^2 = <Hess S v,Hess S v>; in the 1-D slice that is H_vv^2.")


# ---------------------------------------------------------------- [3] the rigorous full-space floor
P("\n[3] full space, NO approximation:  |grad G_v|^2 = <Hess S v, Hess S v> >= (v^T Hess S v)^2 = H_vv^2")
EH_meas = 6.331454          # w7/p02_dyndyn.log line 13 (L=2 SU(2) beta=2.5, combined)
VarG_meas = 5.447690        # w7/p02_dyndyn.log line 15
VarP_meas = 3.61e-3         # w7/o09_phicurv.log line 43 (a=0.75 bin; CONDITIONED, single plaquette)
cs = INV_LAM * EH_meas ** 2
P("    E|grad G_v|^2 >= E[H_vv^2] >= (E H_vv)^2 = %.4f^2 = %.3f" % (EH_meas, EH_meas ** 2))
P("    => (12/N) E|grad G_v|^2 >= %.3f   vs   E[H_vv] = %.4f   ->  ratio %.1f x" %
  (cs, EH_meas, cs / EH_meas))
P("    => line (a): (12/N)E|grad G_v|^2 <= E[H_vv] is FALSE by >= %.0f x  (not a 1-Cov gap)" % (cs / EH_meas))
P("    lambda needed for the BL bound to fall to E[H_vv]:  lam >= E|grad G|^2/E[H_vv] >= E[H_vv] = %.3f"
  % EH_meas)
P("    N/12 = %.4f  -> short by %.1f x   => line (b): BL at N/12 is far too weak." % (LAM, EH_meas / LAM))

# ---------------------------------------------------------------- [4] compare with the measurements
P("\n[4] bound vs measured")
for th in (0.0, np.pi / 2, np.pi):
    m = rows[th]
    P("    th=%.3f: BL(G_v) = %9.3f  vs measured Var(G_v)=%.4f (L=2)  ratio %8.1f x ;"
      "  BL(Phi') = %9.3f vs Var(Phi')=%.2e (cond. a=0.75)  ratio %.1e x"
      % (th, INV_LAM * m["EH2"], VarG_meas, INV_LAM * m["EH2"] / VarG_meas,
         INV_LAM * m["EP22"], VarP_meas, INV_LAM * m["EP22"] / VarP_meas))
P("    caveat: BL bounds the UNconditioned variance; Var(Phi')=3.6e-3 is a CONDITIONED (a=0.75) bin")
P("    value and Var(G_v)=5.45 is the constrained xi=0 slice value -> the comparison is upper-side only.")

# ---------------------------------------------------------------- [5] verdict
P("\n[5] verdict")
P("    (a) NOT satisfied: (12/N)E|grad G_v|^2 >= 240.6 > E[H_vv] = 6.331 (by >=38x).")
P("    (b) YES: BL at lam=N/12 is too weak by 1-3 orders of magnitude; needed lam ~ E[H_vv] ~ 6.3")
P("        (for Phi': needed lam ~ E|grad Phi'|^2/Var(Phi') ~ 1e3).")
P("    (c) YES, and it is the decisive one: Hess(S_W+Phi) is NOT >= (N/12)I at beta=2.5")
P("        (Savvidy: min H_vv=-1.893 pointwise; full Hessian min~-2.8). BL's hypothesis fails.")
P("    => the N/12 convexity CANNOT be converted into a variance bound for the Wilson measure;")
P("       what it proves is the beta=0 (pure Haar) statement, where G_v == 0 identically.")
P("=" * 96)

with open("w7/q02_bl_phi.log", "w", encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
