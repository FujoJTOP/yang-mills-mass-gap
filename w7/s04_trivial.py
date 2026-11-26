"""
S04-TRIVIAL -- is (★) trivially true?

(★)  Var_{mu(.|U')}(<P^T grad S_j, v>) <= E_{mu(.|U')}<P^T H_j P v, v> + W_j(v) - delta_W
 <=> G''_j = E[H_j] - Var_j + W_j >= delta_W > 0

Model used for the exact numbers: single plaquette SU(2), flat log-chart, slice a=0
(xi_3 = 0), v = e_3, boundary holonomy theta.  (H1A-EXACT / O05-FORCEVAR convention.)

    xi in R^3, rho=|xi|, u = xi_3/rho
    S(xi)   = beta (1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) u)
    J        = sin^2(rho/2)/rho^2            (Haar Jacobian),  Phi = -log J
    slice measure  w(rho) drho = 2 pi sin^2(rho/2)/rho * exp(-beta(1-cos(rho/2)cos(th/2))) drho

On the slice (u=0), exactly:
    G   = d_3 S       = beta sin(th/2) sin(rho/2)/rho
    H   = d_3^2 S     = beta cos(th/2) sin(rho/2)/(2 rho)          (theta=pi => H == 0)
    Phi'= 0           (Phi radial) => Var(Phi')=0, X=-2Cov(G,Phi')=0
    Phi''= 2/rho^2 - cot(rho/2)/rho                                (rho->0: 1/6 = N/12, N=2)

Decomposition (G4-GENERAL-BOUND (1.1), convention R1: H = Hess S_W):
    Delta_dyn = E[H] - Var(G) ,  W_geo = E[Phi''] - Var(Phi') ,  X = -2Cov(G,Phi')
    G''       = Delta_dyn + W_geo + X = E[H] - Var_j + W_j ,  W_j := W_geo + X

Convention R2 (H_j = Hess V = S_W''+Phi''):  E[H_j] = E[H]+E[Phi''], Var_j = Var(G+Phi'),
    W_j = 0, and G'' = E[H_j] - Var_j  (identity Var = E[V''] - Var iff IBP holds).

Pure numpy; 1D Gauss-Legendre quadrature, deterministic.
"""
import numpy as np

OUT = []


def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s)
    print(s, flush=True)


# ---------------------------------------------------------------- quadrature ------
RHO_MAX = 2.0 * np.pi


def gl(n=4000, a=0.0, b=RHO_MAX):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (a + b), 0.5 * (b - a) * w


XQ, WQ = gl()


def moments(beta, th):
    """Exact 1D moments on the a=0 slice. Returns dict."""
    c, s = np.cos(th / 2.0), np.sin(th / 2.0)
    rho = XQ
    sr = np.sin(rho / 2.0)
    # slice measure (unnormalised), regularised endpoints
    w = 2.0 * np.pi * sr ** 2 / rho * np.exp(-beta * (1.0 - c * sr))
    G = beta * s * sr / rho
    H = beta * c * sr / (2.0 * rho)
    Phi2 = 2.0 / rho ** 2 - np.cos(rho / 2.0) / sr / rho       # 2/rho^2 - cot(rho/2)/rho
    Z = np.sum(WQ * w)
    E = lambda f: np.sum(WQ * w * f) / Z
    EH, EG, EG2, EP2 = E(H), E(G), E(G ** 2), E(Phi2)
    VarG = EG2 - EG ** 2
    VarPhi = 0.0                      # Phi' == 0 on the slice (exact)
    Xcov = 0.0                        # Cov(G,Phi') = 0 (exact)
    W_geo = EP2 - VarPhi
    Ddyn = EH - VarG
    Gpp = Ddyn + W_geo + Xcov
    return dict(beta=beta, th=th, Z=Z, EH=EH, EG=EG, VarG=VarG, EP2=EP2,
                VarPhi=VarPhi, W_geo=W_geo, X=Xcov, Ddyn=Ddyn, Gpp=Gpp,
                Var_j=VarG + VarPhi,
                EH_R2=EH + EP2, W_j_R1=W_geo + Xcov)


def header(t):
    say("")
    say("=" * 96)
    say(t)
    say("=" * 96)


# ---------------------------------------------------------------- sec 0 -----------
def sec0():
    header("§0  FREE-FIELD / ACID CHECKS and quadrature sanity")
    # W_0 = E[Phi''] at beta=0 must be theta-independent (Haar) and = 1/(2 Cin(2pi))
    m0 = moments(0.0, np.pi)
    say(f"  beta=0, theta=pi : E[Phi''] = {m0['EP2']:.12f}   (target 1/(2 Cin(2pi)) = 0.205115297123)")
    for th in (0.0, np.pi / 2, np.pi):
        m = moments(0.0, th)
        say(f"  beta=0, theta={th:.4f} : E[Phi''] = {m['EP2']:.12f} , E[H] = {m['EH']:.3e} , Var(G) = {m['VarG']:.3e}")
    say("  => beta=0 is theta-INDEPENDENT (Haar); E[H]=0 and Var(G)=0 exactly (G,H ~ beta).")
    # theta=pi for all beta: E[H] identically 0, Var(G) = beta^2 V0
    say("")
    say("  theta=pi family (Haar slice): Var(G)=beta^2*V0 with V0 beta-independent?")
    V0s = []
    for beta in (0.5, 1.0, 2.5, 6.0, 10.0):
        m = moments(beta, np.pi)
        V0s.append(m["VarG"] / beta ** 2)
        say(f"    beta={beta:5.2f}: E[H]={m['EH']:+.2e}  Var(G)/beta^2={m['VarG']/beta**2:.12f}  E[Phi'']={m['EP2']:.9f}")
    say(f"  => V0 = {np.mean(V0s):.12f} (spread {np.ptp(V0s):.2e}) ; target 0.011059293772974")


# ---------------------------------------------------------------- sec 1 -----------
def sec1():
    header("§1  THE QUESTION: beta=0 (pure Haar).  E[H_j], Var_j, W_j, G''_j")
    say("  convention R1 (H_j = Hess S_W;  W_j = W_geo + X)  and  R2 (H_j = Hess V; W_j = 0)")
    say("")
    say(f"  {'theta':>7} {'E[H_j]R1':>12} {'Var_j':>12} {'W_j(R1)=W_geo':>15} {'G\"_j':>12} | "
        f"{'E[H_j]R2':>12} {'W_j(R2)':>8}")
    for th in (0.0, np.pi / 2, np.pi):
        m = moments(0.0, th)
        say(f"  {th:7.4f} {m['EH']:12.6e} {m['Var_j']:12.6e} {m['W_j_R1']:15.9f} {m['Gpp']:12.9f} | "
            f"{m['EH_R2']:12.9f} {0.0:8.1f}")
    say("")
    say("  [T] beta=0 is exact, not numerical:  S_W = beta*S_W^0  =>  G,H == 0 identically,")
    say("      and Phi radial => Phi' == 0 on the slice => Var(Phi')=0, X=0.")
    say("      So  G''(beta=0) = W_geo = E[Phi''] = 0.205115297 = 1/(2Cin(2pi))  (theta-free).")
    say("  => the 0.205 is W_j / E[Phi''], NOT E[H_j]:  E[H_j](beta=0) = 0.")
    say("  => delta_W at beta=0 can be as large as 0.2051  (slack = W_geo, 2.05x the hoped 0.1).")


# ---------------------------------------------------------------- sec 2 -----------
def sec2():
    header("§2  beta>0: does W_j go negative?  Can it cancel E[H_j]?  (exact solvable family)")
    say(f"  {'beta':>8} {'theta':>7} {'E[H]':>10} {'Var_j':>10} {'W_geo':>10} {'X':>6} {'W_j':>10} "
        f"{'Ddyn':>10} {'G\"':>11} {'G\"-0.1':>10}")
    for th in (np.pi, 0.0):
        for beta in (0.0, 0.5, 1.0, 2.5, 3.0, 3.5, 4.0, 4.30660785, 5.0, 6.0, 10.0, 20.0):
            m = moments(beta, th)
            say(f"  {beta:8.4f} {th:7.4f} {m['EH']:+10.4f} {m['Var_j']:10.6f} {m['W_geo']:10.6f} "
                f"{m['X']:+6.1f} {m['W_j_R1']:+10.6f} {m['Ddyn']:+10.4f} {m['Gpp']:+11.6f} {m['Gpp']-0.1:+10.4f}")
        say("")

    # exact closed forms at theta=pi
    beta = np.linspace(0.0, 8.0, 16001)
    g = np.array([moments(b, np.pi)["Gpp"] for b in beta[::40]])
    bb = beta[::40]
    idx = np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]
    for i in idx:
        lo, hi = bb[i], bb[i + 1]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if moments(lo, np.pi)["Gpp"] * moments(mid, np.pi)["Gpp"] <= 0:
                hi = mid
            else:
                lo = mid
        say(f"  root  G''(beta*)=0  at beta* = {0.5*(lo+hi):.9f}   (archive: 4.30660785)")
    # delta_W = 0.1 window edge: G''(beta) = 0.1  (own bracket search)
    h = np.array([moments(b, np.pi)["Gpp"] - 0.1 for b in bb])
    for i in np.where(np.sign(h[:-1]) != np.sign(h[1:]))[0]:
        lo2, hi2 = bb[i], bb[i + 1]
        for _ in range(60):
            mid = 0.5 * (lo2 + hi2)
            if (moments(lo2, np.pi)["Gpp"] - 0.1) * (moments(mid, np.pi)["Gpp"] - 0.1) <= 0:
                hi2 = mid
            else:
                lo2 = mid
        say(f"  root  G''(beta)=0.1 at beta  = {0.5*(lo2+hi2):.9f}   (delta_W=0.1 window edge)")
    W0, V0 = moments(0.0, np.pi)["W_geo"], moments(1.0, np.pi)["VarG"]
    say(f"  closed form: beta*(delta_W) = sqrt((W0-delta_W)/V0) = {np.sqrt((W0-0.1)/V0):.9f}")


# ---------------------------------------------------------------- sec 3 -----------
def sec3():
    header("§3  Non-degenerate slices (cos(theta/2) != 0): does (★) survive large beta?")
    say("  O05-FORCEVAR claim: m=dlogE[H]/dlogbeta -> 1, k=dlogVar/dlogbeta -> 0, so RHS wins at large beta.")
    say(f"  {'theta':>7} {'beta':>7} {'E[H]':>11} {'Var_j':>11} {'W_geo':>10} {'Ddyn':>11} {'G\"':>12}")
    for th in (0.0, np.pi / 3, np.pi / 2):
        for beta in (2.5, 6.0, 20.0, 50.0, 200.0):
            m = moments(beta, th)
            say(f"  {th:7.4f} {beta:7.1f} {m['EH']:11.5f} {m['Var_j']:11.6f} {m['W_geo']:10.6f} "
                f"{m['Ddyn']:11.5f} {m['Gpp']:12.6f}")
        say("")
    # local powers
    say("  local powers (beta 50 -> 200):")
    for th in (0.0, np.pi / 2):
        a, b = moments(50.0, th), moments(200.0, th)
        k = np.log(max(b["Var_j"], 1e-300) / max(a["Var_j"], 1e-300)) / np.log(4.0)
        mm = np.log(max(b["EH"], 1e-300) / max(a["EH"], 1e-300)) / np.log(4.0)
        say(f"    theta={th:.4f}: k=Var power={k:.4f} , m=E[H] power={mm:.4f} , "
            f"E[H]/beta={a['EH']/50:.6f}/{b['EH']/200:.6f}")


# ---------------------------------------------------------------- sec 4 -----------
def sec4():
    header("§4  Can W_j (not Var_j) cancel E[H_j]?  Cauchy-Schwarz ceiling on X = -2Cov(G,Phi')")
    say("  On the a=0 SU(2) slice X == 0 EXACTLY (Phi radial => Phi' == 0).  So W_j is pinned at")
    say("  W_geo = 0.2051 (theta=pi, beta-free) and CANNOT cancel anything there: the cancellation is")
    say("  done 100% by Var_j = beta^2 V0.")
    say("")
    say("  Off the slice X != 0.  Universal ceiling: |X| <= 2 sqrt(VarG * VarPhi')  (Cauchy-Schwarz).")
    say("  Measured inputs (archive, SU(3) L=4/6 pool, beta=5.7 / 2.5):")
    say("    VarG ~ 2.0 ,  VarPhi' ~ 0.02  =>  2 sqrt(0.04) = 0.400   vs   W_geo ~ 0.32")
    say("    => W_j >= 0.32 - 0.40 = -0.08  < 0   (allowed, and 0.08 ~ the whole N/12=0.167 floor/2)")
    say("    measured W_full = W_geo + X ~ 0.54 > 0.32  =>  measured X ~ +0.22 (POSITIVE at that point)")
    say("  => W_j CAN be negative enough to eat E[H_j] on general slices (there Ddyn = E[H]-VarG ~ 0,")
    say("     equipartition), but on the exactly solvable slice it is Var_j, not W_j, that does the eating.")
    say("     Either way (★) is a real inequality with a real failure window -- not a tautology.")


# ---------------------------------------------------------------- sec 5 -----------
def sec5():
    header("§5  VERDICT INPUTS")
    m0 = moments(0.0, np.pi)
    say(f"  beta=0 :  E[H_j](R1)=0 ,  Var_j=0 ,  W_j=W_geo={m0['W_geo']:.9f} ,  G''={m0['Gpp']:.9f}  [T]+[V]")
    say(f"           (R2) E[H_j]={m0['EH_R2']:.9f}=E[Phi''] , W_j=0 ,  G'' same")
    say("  beta>0 , theta=pi :  G''(beta) = W_geo - beta^2 V0 = 0.205115297 - 0.011059294 beta^2")
    say("           G''=0 at beta*=4.30661 ; G''=0.1 at beta=3.0776 ;  G''<0 for beta>beta*  [C]")
    say("  => (★) holds at beta=0 with slack 0.205 (delta_W=0.1 free), is NOT trivial, and FAILS")
    say("     for beta>beta* on the a0=0, theta=pi, v=e3 family (= V3/P18 [NG] witness, reproduced).")


if __name__ == "__main__":
    sec0()
    sec1()
    sec2()
    sec3()
    sec4()
    sec5()
    open(__file__.replace(".py", ".log"), "w", encoding="utf-8").write("\n".join(OUT) + "\n")
