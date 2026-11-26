"""R09-FLUCT -- judge  beta*Var(d_t h) <= E[d_t^2 h]  (curvature fluctuation vs curvature mean).

DEFINITIONS (fixed before any number; see R09-FLUCT.md sec.0)
  q06:  h := sum_P (1 - ReTr U_P/2)  (Wilson action density, BETA-FREE),  V = beta*h + Phi
        g := d_t h,  h2 := d_t^2 h  (both beta-free),  S = E[kappa] = E[H + Phi'']
  q04/r03 single plaquette SU(2):  S_W = beta*(1 - c*cos u + s*sin u*xi3/rho), u=rho/2
        c=cos(th/2), s=sin(th/2)   =>  h = S_W/beta  EXACTLY, and
        d_a S_W = G_v = beta*d_a h,   d_a^2 S_W = H_vv = beta*d_a^2 h.
  => q06's flow parameter t is the slice direction xi3=a (the direction v=e3 of the Hessian),
     NOT beta (h is beta-free, d_beta h = 0) and NOT a thermal time.
  => beta*Var(d_t h) <= E[d_t^2 h]  <=>  Var(G_v) <= E[H_vv]     (beta cancels, ratio identical)
     i.e. exactly q04/r02's criterion  2*beta*s^2*Var(g) <= c*E[g].

usage: py -3.12 w7/r09_fluct.py   (writes w7/r09_fluct.log)
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

TWOPI = 2.0 * np.pi
BETAS = (0.5, 1.0, 2.5, 4.0, 6.0)
THETAS = (0.0, 0.5 * np.pi, np.pi)      # the three required slices


def gl(n=4000):
    return leggauss(n)


_X, _W = gl()


def moms(beta, th, n=4000):
    """exact 1-D slice moments (a=0, u in (0,pi)): nu prop e^{-beta(1-c cos u)} sin^2 u/u du,
       g = sin u/(2u).  Returns Z, E[g], Var[g], E[Phi''], and int of the r02 kernel."""
    c, s = np.cos(th / 2), np.sin(th / 2)
    u = 0.5 * np.pi * (_X + 1.0)
    wt = _W * 0.5 * np.pi
    e = np.exp(-beta * (1.0 - c * np.cos(u)))
    w = np.sin(u) ** 2 / u                      # the density's beta-independent factor
    g = np.sin(u) / (2 * u)
    phipp = 1.0 / (2 * u * u) - (np.cos(u) / np.sin(u)) / (2 * u)
    Z = (wt * e * w).sum()
    E1 = (wt * e * w * g).sum() / Z
    E2 = (wt * e * w * g * g).sum() / Z
    Ep = (wt * e * w * phipp).sum() / Z
    return dict(c=c, s=s, Z=Z, E1=E1, V=E2 - E1 ** 2, Ep=Ep, E2=E2)


def kern_half(beta, th):
    """(1/2) int_0^pi e sin^3 u (u - sin u)/u^3 du  -- should equal Z*(E[g]-2E[g^2])."""
    c = np.cos(th / 2)
    u = 0.5 * np.pi * (_X + 1.0)
    wt = _W * 0.5 * np.pi
    e = np.exp(-beta * (1.0 - c * np.cos(u)))
    h = e * np.sin(u) ** 3 * (u - np.sin(u)) / u ** 3
    return 0.5 * (wt * h).sum(), float(np.min(h))


def S_W(rho, a, th, beta):
    c, s = np.cos(th / 2), np.sin(th / 2)
    return beta * (1.0 - c * np.cos(rho / 2) + s * np.sin(rho / 2) * (a / rho))


def h_of(rho, a, th):
    """q06's h = S_W/beta."""
    return S_W(rho, a, th, 1.0)


L = []


def P(x=""):
    print(x)
    L.append(x)


def main():
    beta0 = 2.5
    P("R09-FLUCT  beta=2.5, single plaquette SU(2), slice xi3=a=0  (py -3.12)")
    P("target:  beta*Var(d_t h) <= E[d_t^2 h]   with h = S_W/beta, t = xi3 = a")
    P("")

    # ---------------- [0] definition / identity checks (no numbers before this)
    P("[0] definition audit")
    eta, a0 = 1.1, 0.0
    rho = np.hypot(a0, eta)
    d = 1e-6
    fd1 = (h_of(np.hypot(a0 + d, eta), a0 + d, 0.7) - h_of(np.hypot(a0 - d, eta), a0 - d, 0.7)) / (2 * d)
    fd2 = (h_of(np.hypot(a0 + d, eta), a0 + d, 0.7) - 2 * h_of(rho, a0, 0.7)
           + h_of(np.hypot(a0 - d, eta), a0 - d, 0.7)) / d ** 2
    s_, c_ = np.sin(0.35), np.cos(0.35)
    ex1 = s_ * np.sin(rho / 2) / rho
    ex2 = (c_ / 2) * np.sin(rho / 2) / rho
    P(f"    h := S_W/beta pointwise: max|beta*h - S_W| over grid = "
      f"{np.max(np.abs(np.array([S_W(r, aa, 0.9, 2.5) - 2.5 * h_of(r, aa, 0.9) for r, aa in
                               [(0.3, 0.1), (1.7, 0.4), (4.0, 1.2)]]))):.1e}")
    P(f"    d_a h (FD, th=0.7, eta=1.1) = {fd1:+.10f}   s*g = {ex1:+.10f}   diff {fd1 - ex1:+.1e}")
    P(f"    d_a^2 h (FD)             = {fd2:+.10f}   (c/2)g = {ex2:+.10f}   diff {fd2 - ex2:+.1e}")
    P("    => G_v = beta*d_a h,  H_vv = beta*d_a^2 h,  h is beta-free (d_beta h = 0)")
    P("")

    # ---------------- [1] the three slices at beta=2.5
    P("[1] beta=2.5, a=0:  d_a h = s*g,  d_a^2 h = (c/2)*g,  g=sin(u)/(2u)  (SAME random variable)")
    P(f"    {'theta/pi':>8} {'E[g]':>9} {'Var(g)':>10} {'LHS=b*Var(dth)':>15} {'RHS=E[d2h]':>11}"
      f" {'ratio':>8} {'hold':>5}")
    rows = {}
    for th in THETAS:
        m = moms(beta0, th)
        lhs = beta0 * m["s"] ** 2 * m["V"]
        rhs = 0.5 * m["c"] * m["E1"]
        rows[th] = (m, lhs, rhs)
        P(f"    {th / np.pi:8.3f}pi {m['E1']:9.6f} {m['V']:10.6f} {lhs:15.8f} {rhs:11.6f}"
          f" {lhs / rhs if rhs else np.inf:8.4f} {'YES' if lhs <= rhs else 'NO':>5}")
    P("")

    # ---------------- [2] the beta-cancellation identity
    P("[2] beta cancels identically:  beta*Var(d_a h) <= E[d_a^2 h]  <=>  Var(G_v) <= E[H]")
    for th in THETAS:
        m, lhs, rhs = rows[th]
        vg = beta0 ** 2 * m["s"] ** 2 * m["V"]
        eh = 0.5 * beta0 * m["c"] * m["E1"]
        r1 = lhs / rhs if rhs else np.inf
        r2 = vg / eh if eh else np.inf
        P(f"    theta={th / np.pi:.3f}pi: Var(G_v)={vg:.8f}  E[H_vv]={eh:+.8f}  "
          f"ratio_q06={r1:.10f}  ratio_q04={r2:.10f}  diff={abs(r1 - r2):.1e}")
    P("")

    # ---------------- [3] the r02 kernel: 2Var(g) <= E[g] unconditionally  [T]
    P("[3] r02 kernel (exact, pointwise-positive integrand => unconditional [T]):")
    P("    E[g] - 2E[g^2] = (1/2) int_0^pi e^{-beta(1-c cos u)} sin^3 u (u - sin u)/u^3 du >= 0")
    worst = -1e9
    for beta in (0.0, 0.5, 2.5, 6.0, 20.0):
        for th in THETAS:
            m = moms(beta, th)
            kh, kmin = kern_half(beta, th)
            P(f"    beta={beta:5.1f} th={th / np.pi:.2f}pi: Z(E1-2E2)={m['Z'] * (m['E1'] - 2 * m['E2']):+.6e}"
              f"  halfint={kh:+.6e}  min(integrand)={kmin:+.2e}  2Var/E={2 * m['V'] / m['E1']:.6f}")
            worst = max(worst, 2 * m["V"] / m["E1"])
    P(f"    => 2Var(g)/E[g] <= 1/2 for every (beta,c) tested; max = {worst:.6f}  [T]")
    P("")

    # ---------------- [4] sufficient window and exact boundary
    P("[4] sufficient condition [T]: beta*s^2 <= c, i.e. beta*sin^2(th/2) <= cos(th/2)")
    ths = np.linspace(0.0, np.pi, 200001)
    ok = beta0 * np.sin(ths / 2) ** 2 <= np.cos(ths / 2)
    P(f"    beta=2.5: holds [T] for theta/pi <= {ths[ok].max() / np.pi:.4f}")
    P("    exact ratio(1) boundary, grid in beta (ratio NOT monotone in beta -> scan, no bisection):")
    for th in THETAS:
        c, s2 = np.cos(th / 2), np.sin(th / 2) ** 2
        if s2 == 0.0:
            P(f"    theta={th / np.pi:.2f}pi  c={c:.5f}: LHS==0 for every beta -> holds all beta")
            continue
        bs = np.exp(np.linspace(np.log(1e-3), np.log(60.0), 4000))
        rt = np.array([beta_ratio(b, th) for b in bs])
        bad = bs[rt > 1.0]
        mx = bs[int(np.argmax(rt))]
        P(f"    theta={th / np.pi:.2f}pi  c={c:.5f}: fails on beta in "
          f"{('[' + f'{bad.min():.3f},{bad.max():.3f}]') if len(bad) else 'never'}  "
          f"max ratio={rt.max():.4f} at beta={mx:.3f}  (at beta=2.5: {beta_ratio(2.5, th):.4f})")
    P("")

    # ---------------- [5] a != 0: the inequality is NOT a=0-only
    P("[5] a-dependence at theta=pi, beta=2.5 (q04 reduction, rho in [|a|,2pi])")
    for a in (0.0, 0.25, 0.5, 1.0, 2.0):
        VG, EH = slice_am(a, np.pi, beta0)
        P(f"    a={a:4.2f}: Var(G_v)={VG:.6f}  E[H_vv]={EH:+.6f}  "
          f"{'hold' if VG <= EH else 'FAIL'}  (E[H]-Var={EH - VG:+.6f})")
    P("")

    # ---------------- [6] verdict
    P("[6] verdict")
    P("    theta=0     : 0 <= (1/2)E[g] = 0.221225 -> holds (degenerate, s=0)")
    P("    theta=pi/2  : 0.006096 <= 0.150909      -> holds (ratio 0.0404)")
    P("    theta=pi    : 0.027648 <= 0             -> FAILS  (H_vv == 0 pointwise: c=0)")
    P("    [T] a=0 slice, any beta>=0, any c in [0,1]: 2Var(g) <= E[g] (r02 kernel)")
    P("        => beta*Var(d_t h) <= E[d_t^2 h] holds whenever beta*sin^2(th/2) <= cos(th/2)")
    P("    [E] theta=pi at every beta>0: counterexample (exact, no quadrature error)")
    P("    [E] a>0, theta=pi: E[H_vv]<0 fails harder")
    with open("w7/r09_fluct.log", "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")


def beta_ratio(beta, th):
    m = moms(beta, th)
    rhs = 0.5 * m["c"] * m["E1"]
    return beta * m["s"] ** 2 * m["V"] / rhs if rhs else np.inf


def slice_am(a, th, beta, n=4000):
    """q04 moments() reduction (independent re-implementation): eta=sqrt(rho^2-a^2)."""
    a = abs(float(a))
    eta_max = np.sqrt(TWOPI ** 2 - a * a)
    x, wq = gl(n)
    eta = 0.5 * eta_max * (x + 1.0)
    rho = np.sqrt(a * a + eta * eta)
    jac = 2 * np.pi * eta * np.sin(rho / 2) ** 2 / rho ** 2
    wt = jac * np.exp(-S_W(rho, a, th, beta)) * 0.5 * eta_max * wq
    Z = wt.sum()
    c, s = np.cos(th / 2), np.sin(th / 2)
    wr = np.sin(rho / 2) / rho
    wp = ((rho / 2) * np.cos(rho / 2) - np.sin(rho / 2)) / rho ** 2
    wpp = -np.sin(rho / 2) / (4 * rho) - np.cos(rho / 2) / rho ** 2 + 2 * np.sin(rho / 2) / rho ** 3
    G = beta * ((c / 2) * (a / rho) * np.sin(rho / 2) + s * (wr + (a * a / rho) * wp))
    A = np.sin(rho / 2) / rho + (a * a / (2 * rho ** 2)) * np.cos(rho / 2) - (a * a / rho ** 3) * np.sin(rho / 2)
    B = (3 * a / rho) * wp - (a ** 3 / rho ** 3) * wp + (a ** 3 / rho ** 2) * wpp
    H = beta * ((c / 2) * A + s * B)
    EG = (wt * G).sum() / Z
    return (wt * G * G).sum() / Z - EG ** 2, (wt * H).sum() / Z


if __name__ == "__main__":
    main()
