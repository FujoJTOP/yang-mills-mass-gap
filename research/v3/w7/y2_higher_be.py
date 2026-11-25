#!/usr/bin/env python3
r"""
y2_higher_be.py -- Y2: test the HIGHER-ORDER / NONLOCAL Bakry-Emery proposal against the Savvidy
soft mode.  Standard BE/LSI dies because G'' (the marginal curvature = the one-sided BE constant,
H6-B) has ess inf < 0 on the Savvidy tail.  Question: can a criterion using HIGHER covariant
derivatives (iterated Gamma_k / k-convexity / regularized curvature V''+a V'''') or a NONLOCAL
functional restore a uniform log-concavity => Poincare bound?

Tests:
 (A) FIXED-JET DOUBLE WELL [T]+[V].  V_w(y)=(w^2/2)((y/w)^2-1)^2 : barrier curvature V''(0)=-2,
     well curvature V''(+-w)=+4 are INDEPENDENT of w; V''''(0)=12/w^2 -> 0; barrier height w^2/2
     -> infinity.  => the FULL jet at every critical point converges to a fixed jet while lambda_1 -> 0.
     Hence NO criterion that is a function of the critical-point jets can give a uniform floor.
 (B) 4-CONVEX / REGULARIZED CURVATURE [T]+[V].  V_g=(g/2)(x^2-1)^2 is 4-convex (V''''=6g>0) and
     the regularized curvature kappa_a=V''+a V'''' is +2g>0 at the barrier for a=1/3 (and +4g for a=1),
     YET lambda_1 -> 0 (Kramers).  => "higher-order log-concavity" is satisfied but gives no gap.
 (C) FREE FIELD [T]+[V].  V=1/2 k^2 x^2 : every jet of order >=3 is 0; any higher-order criterion
     returns 0; lambda_1=k^2 -> 0.  ACID PASS.
 (D) EXACT ONE-PLAQUETTE SU(2) [T]+[V].  G''(e3;th=pi,a0=0)=A-B beta^2, A=0.20511530, B=0.01105929,
     beta*=4.30661.  Compute the exact marginal rho(a), G''(0), G''''(0), lambda_1(nu), Cheeger gate,
     below (beta=2.5) and above (5,10) beta*.  Report whether a higher-order regularization restores a
     positive curvature at the barrier and whether it tracks the actual gap.
 (E) 1D SU(2) EXACT HAAR [T]+[V].  log-concave single well; lambda_1=0.2645, Var=3.7198, chi<=3.78
     (beta=2.5 floor 0.222 => chi<=4.5).  Higher-order is unnecessary; BE already suffices.
 (F) NONLOCAL [T]/[NG].  The only nonlocal functional that controls lambda_1 is the Cheeger/gate
     capacity h(nu) ~ rho(gate) ~ e^{-DeltaV}; it IS H1 (the gate-density bound).

Run: py -3.12 -I y2_higher_be.py .  numpy+scipy, deterministic, small 1-D jobs.  No QED.  Clay 0%.
"""
import numpy as np
from scipy.special import sici
from scipy.linalg import eigh

EULER = 0.5772156649015328606
TWO_PI = 2.0 * np.pi


# ------------------------------------------------------------------ shared 1-D tools
def lam1_of_marginal(x, logrho):
    """Lowest positive eigenvalue of L = d^2 - V' d (V=-log rho) on L^2(nu), P1 FEM. Neumann ends."""
    x = np.asarray(x, float); logrho = np.asarray(logrho, float); n = len(x); h = np.diff(x)
    rho = np.exp(logrho - logrho.max()); rho /= np.trapezoid(rho, x)
    K = np.zeros((n, n)); M = np.zeros((n, n))
    for e in range(n - 1):
        he = h[e]; w = 0.5 * (rho[e] + rho[e + 1]) * he
        ke = (1.0 / he) * np.array([[1.0, -1.0], [-1.0, 1.0]])
        me = (he / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])
        idx = [e, e + 1]
        for ii in range(2):
            for jj in range(2):
                K[idx[ii], idx[jj]] += ke[ii, jj] * w
                M[idx[ii], idx[jj]] += me[ii, jj] * w
    vals = np.sort(eigh(K, M, eigvals_only=True))
    pos = vals[vals > 1e-8]          # exclude the exact kernel 0 and round-off spurious ~1e-11 modes
    return (pos[0] if len(pos) else 0.0)


def cheeger_h(x, rho):
    """Exact 1-D Cheeger constant h = inf_x rho(x)/min(F,1-F).  Rigorous: lambda_1 >= h^2/4.
    (Cheeger constant with h=inf Per/min(mu,mu^c); in 1-D Per((x,inf))=rho(x), no extra factor.)"""
    x = np.asarray(x, float); r = np.asarray(rho, float); r = r / np.trapezoid(r, x)
    F = np.concatenate([[0.0], np.cumsum(0.5 * (r[1:] + r[:-1]) * np.diff(x))]); F = np.clip(F, 0, 1)
    mn = np.minimum(F, 1.0 - F)
    ok = (F > 1e-3) & (F < 1.0 - 1e-3)
    psi = np.where(ok, r / np.maximum(mn, 1e-300), np.inf)
    i = int(np.argmin(psi))
    return psi[i], x[i], F[i]


# ------------------------------------------------------------------ (A) fixed-jet double well
def fixed_jet_family(w, n=3001):
    """V_w(y)=(w^2/2)((y/w)^2-1)^2 on [-1.6w,1.6w].  Barrier V''(0)=-2, well V''(+-w)=+4 fixed
    in w; V''''=12/w^2 -> 0; barrier height w^2/2 -> infinity."""
    y = np.linspace(-1.6 * w, 1.6 * w, n)
    V = (w ** 2 / 2.0) * ((y / w) ** 2 - 1.0) ** 2
    logrho = -V; logrho -= logrho.max()
    Vpp = -2.0 + 6.0 * (y / w) ** 2          # exact V''
    V4 = np.full_like(y, 12.0 / w ** 2)       # exact V''''
    barriers = V[np.argmin(np.abs(y))] - V.min()   # V(0)-V(well) = w^2/2
    lam1 = lam1_of_marginal(y, logrho)
    return lam1, barriers, Vpp[np.argmin(np.abs(y))], Vpp[np.argmin(np.abs(np.abs(y) - w))], V4[0]


# ------------------------------------------------------------------ (B) 4-convex g-family
def g_family(g, n=3001):
    x = np.linspace(-2.6, 2.6, n)
    V = 0.5 * g * (x ** 2 - 1.0) ** 2
    logrho = -V; logrho -= logrho.max()
    Vpp = g * (3.0 * x ** 2 - 1.0)            # V''
    V4 = np.full_like(x, 6.0 * g)             # V''''
    lam1 = lam1_of_marginal(x, logrho)
    return lam1, Vpp.min(), Vpp.max(), V4[0]


# ------------------------------------------------------------------ (D) exact one-plaquette SU(2)
XG, WG = np.polynomial.legendre.leggauss(3000)


def phi0(r):
    r = np.asarray(r, float); sm = r < 1e-6; rs = np.where(sm, 1.0, r)
    val = 2.0 * np.log(rs / 2.0) - 2.0 * np.log(np.sin(rs / 2.0))
    return np.where(sm, r * r / 12.0, val)


def rho_v_1x1(a0, b, t):
    """Exact one-link SU(2) marginal in the e3 direction, boundary holonomy theta=t (S2/G4/V3)."""
    Rm = np.sqrt(max(4 * np.pi ** 2 - a0 * a0, 1e-16)); rr = 0.5 * Rm * (XG + 1.0); wrr = 0.5 * Rm * WG
    rho = np.sqrt(rr * rr + a0 * a0); u = a0 / rho
    V = b * (1.0 - np.cos(rho / 2.0) * np.cos(t / 2.0)) + phi0(rho) + b * np.sin(rho / 2.0) * np.sin(t / 2.0) * u
    return 2.0 * np.pi * np.sum(np.exp(-V) * rr * wrr)


def one_plaquette(beta, th, n=2401):
    a = np.linspace(-TWO_PI + 1e-3, TWO_PI - 1e-3, n)
    rho = np.array([rho_v_1x1(x, beta, th) for x in a])
    rho = np.maximum(rho, 1e-300)
    lr = np.log(rho); lr -= lr.max()
    # derivatives of G''=-d^2 log rho/da^2, G''''=-d^4 log rho/da^4 (finite diff, refined near a=0)
    da = a[1] - a[0]
    G2 = -np.gradient(np.gradient(lr, a), a)
    G4 = -np.gradient(np.gradient(G2, a), a)
    lam1 = lam1_of_marginal(a, lr)
    hh, xg, Fg = cheeger_h(a, rho)
    # robust central stencil near a=0
    i0 = int(np.argmin(np.abs(a)))
    def d2(i, f, h): return (-f[i - 2] + 16 * f[i - 1] - 30 * f[i] + 16 * f[i + 1] - f[i + 2]) / (12 * h * h)
    def d4(i, f, h):
        return (f[i - 2] - 4 * f[i - 1] + 6 * f[i] - 4 * f[i + 1] + f[i + 2]) / (h ** 4)
    Gpp0 = -d2(i0, lr, da)
    G4p0 = -d4(i0, lr, da)
    a_star = a[int(np.argmax(rho))]
    return dict(lam1=lam1, Gpp0=Gpp0, G4p0=G4p0, h=hh, gate=xg, a_star=a_star,
                G2_curve=(a, G2), rho=(a, rho))


# ------------------------------------------------------------------ (E) 1D SU(2) exact Haar
def Cin(x):
    x = np.asarray(x, float); _, ci = sici(x)
    out = EULER + np.log(np.where(x > 0, x, 1.0)) - ci
    return np.where(x > 0, out, 0.0)


def haar_marginal():
    z = np.linspace(-TWO_PI, TWO_PI, 1601)[1:-1]
    rho = np.pi * (Cin(TWO_PI) - Cin(np.minimum(np.abs(z), TWO_PI)))
    lr = np.log(np.maximum(rho, 1e-300)); lr -= lr.max()
    rho = np.exp(lr); rho /= np.trapezoid(rho, z)
    var = np.trapezoid((z - np.trapezoid(z * rho, z)) ** 2 * rho, z)
    lam1 = lam1_of_marginal(z, lr)
    return var, lam1, 1.0 / var


# ------------------------------------------------------------------ main
def main():
    print("=" * 90)
    print("(A) FIXED-JET DOUBLE WELL V_w=(w^2/2)((y/w)^2-1)^2 : jets at critical points fixed in w,")
    print("    barrier height w^2/2 -> inf, lambda_1 -> 0  [T]+[V]")
    print("      w    V''(barrier)   V''(well)   V''''(barrier)   barrier   lambda_1(FEM)   e^{-barrier}")
    for w in (1.0, 1.5, 2.0, 3.0, 4.0):
        lam1, bar, vpb, vpw, v4 = fixed_jet_family(w)
        print(f"    {w:4.1f}    {vpb:+8.4f}     {vpw:+8.4f}     {v4:11.4e}     {bar:7.3f}   "
              f"{lam1:11.3e}    {np.exp(-bar):9.2e}")
    print("    => V''(0)=-2 and V''(well)=+4 are INDEPENDENT of w; V''''->0. The full critical-point")
    print("       jet converges, yet lambda_1->0. NO criterion f(V'',V''',V'''',...) of the critical-")
    print("       point jets can give a uniform floor: the gap is the NONLOCAL barrier action.")

    print("=" * 90)
    print("(B) 4-CONVEX family V_g=(g/2)(x^2-1)^2 : V''''=6g>0 everywhere, kappa_a=V''+a V'''' at")
    print("    the barrier = (6a-2)g > 0 for a>1/3, yet lambda_1 -> 0 (Kramers)  [T]+[V]")
    print("       g    V''''=6g    V''(barrier)  kappa_1/3(bar)  kappa_1(bar)   lambda_1(FEM)")
    for g in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0):
        lam1, vpmin, vpmax, v4 = g_family(g)
        kap13 = vpmin + (1.0 / 3.0) * v4
        kap1 = vpmin + 1.0 * v4
        print(f"    {g:4.1f}   {v4:9.3f}    {vpmin:+8.3f}     {kap13:+10.3f}     {kap1:+9.3f}    "
              f"{lam1:10.3e}")
    print("    => 'higher-order log-concavity' (V''''>=0) and the regularized curvature are SATISFIED")
    print("       (kappa_1=+4g at the barrier) while the gap collapses. A higher-derivative criterion")
    print("       positive on the double well cannot imply a Poincare bound.  [NG]")

    print("=" * 90)
    print("(C) FREE FIELD V=1/2 k^2 x^2 : all jets of order >=3 vanish; lambda_1=k^2 -> 0  [T]+[V]")
    for k in (1.0, 0.5, 0.1, 0.02):
        x = np.linspace(-6.0 / k, 6.0 / k, 2001)
        lr = -0.5 * k * k * x * x; lr -= lr.max()
        lam1 = lam1_of_marginal(x, lr)
        print(f"    k={k:5.2f}: V''=k^2={k*k:8.4f}  V''''=0  lambda_1(FEM)={lam1:.5f}   (ACID PASS)")
    print("    => every higher-order criterion that is a function of the jet returns 0 on the free")
    print("       field (all higher jets are exactly 0). No false floor is manufactured.")

    print("=" * 90)
    print("(D) EXACT ONE-PLAQUETTE SU(2), theta=pi, G''(0)=A-B beta^2, A=0.20511530 B=0.01105929,")
    print("    beta*=4.30661  [T]+[V]")
    A0, B0 = 0.2051152978, 0.0110592938
    print("     beta   G''(0)     A-B b^2    G''''(0)   kappa1=G''+G''''  lambda_1   h^2/4    gate")
    for beta in (2.5, 5.0, 10.0):
        r = one_plaquette(beta, np.pi)
        pred = A0 - B0 * beta * beta
        kap1 = r['Gpp0'] + r['G4p0']
        print(f"    {beta:5.2f}  {r['Gpp0']:+9.5f}  {pred:+9.5f}  {r['G4p0']:+9.4f}   "
              f"{kap1:+9.4f}        {r['lam1']:9.4f}  {r['h']**2/4:8.4f}  {r['gate']:+6.2f}")
    print("    Below beta* the barrier has G''>0 (log-concave): plain BE works.")
    print("    Above beta* G''(0)<0 but the exact marginal gap stays O(0.4): it is a compact tilted")
    print("    SINGLE well, the G''<0 point is an off-mass shoulder (Q3), not a bottleneck.  Note")
    print("    kappa1=G''+G''''=+0.62 > lambda_1=0.40 at beta=10 (and >lambda_1 everywhere above),")
    print("    so the 'regularized curvature' is NOT a lower bound for lambda_1 on the exact model [NG].")

    print("=" * 90)
    print("(E) 1D SU(2) exact one-link Haar marginal: log-concave single well  [T]+[V]")
    var, lam1, invv = haar_marginal()
    print(f"    G''(0)=1/(2Cin2pi)={1.0/(2*Cin(TWO_PI)):.6f}  Var={var:.4f}  1/Var={invv:.6f}")
    print(f"    lambda_1(nu)={lam1:.6f}  =>  chi <= 1/lambda_1 = {1/lam1:.3f}  (measured Var={var:.3f})")
    print("    ess inf G''>0 (single well): plain BE already suffices; higher-order adds nothing.")
    print("    beta=2.5 periodic k=0 floor (P2/P4): G''->0.222, Var~4.45 => chi<=4.5.")

    print("=" * 90)
    print("(F) NONLOCAL: the only functional controlling lambda_1 is the Cheeger/gate capacity h(nu).")
    print("    On V=(g/2)(x^2-1)^2 the gate density rho(gate)~e^{-DeltaV} and h~rho(gate)/mass:  [T]+[V]")
    print("       g    rho(gate)     Cheeger h    h^2/4      lambda_1     H1?")
    for g in (0.5, 2.0, 8.0, 16.0):
        x = np.linspace(-2.6, 2.6, 2401); V = 0.5 * g * (x ** 2 - 1.0) ** 2
        rho = np.exp(-V); rho /= np.trapezoid(rho, x)
        i0 = int(np.argmin(np.abs(x))); rg = rho[i0]
        hh, _, _ = cheeger_h(x, rho)
        lam1, _, _, _ = g_family(g)
        print(f"    {g:4.1f}  {rg:10.3e}   {hh:10.3e}  {hh**2/4:9.2e}  {lam1:10.3e}   "
              f"gate-bound = H1")
    print("    => h is a NONLOCAL (barrier-action) functional; h>=h0 <=> rho(gate)>=c <=> no interior")
    print("       barrier <=> H1. The nonlocal rescue IS H1. Nonlocal Gamma / fractional LSI still")
    print("       require a pointwise curvature floor (absent). Orderless quotient: half-line optimality")
    print("       of 1-D Cheeger has no analogue on A/G; only the tautology lambda_1>=h^2/4 survives.")

    print("=" * 90)
    print("SUMMARY (no QED):  (A) fixed critical-point jet + lambda_1->0 => no local high-jet criterion;")
    print("  (B) 4-convex / regularized curvature positive but gap 0 => [NG] higher-order log-concavity;")
    print("  (C) free field acid pass; (D) one-plaquette: barrier G''=A-Bb^2<0 past beta*=4.30661;")
    print("  (E) 1D SU(2) finite chi via plain BE; (F) the only working nonlocal functional = Cheeger")
    print("  gate = H1. Uniform bound still [O]=H1. Clay 0%.")


if __name__ == "__main__":
    main()
