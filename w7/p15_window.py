"""
p15_window.py -- P15: does the 3x over-count of E[H] (O15-THETAPI Sec 4.2) propagate into
O21's c and hence into O22's beta0 / window?  Recomputed c, beta0, and the theta-window,
two scenarios (uncorrected / 3x-corrected), sharp (Var<=V_inf) and strict (x6.51) bounds.

Convention: O21's exact single-plaquette SU(2) model, imported directly (no duplication):
    S = beta(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) u) = beta(1 - ReTr(UB)/N), N=2
    slice xi_3=0, v=e_3,  H_vv = beta*cos(th/2)*sin(rho/2)/(2rho),  Phi''=2/rho^2-cot(rho/2)/rho
    c(beta) = E[H_vv+Phi''|0]/(beta||v||^2) = E[Q|0] + E[Phi''|0]/beta
Textbook facts used (all [T]):  c_inf = cos(th/2)/(2N) = Re tr U_p/(4N);
V_inf = tan^2(th/2)/36 (O20 sharp bound); Var <= 6.51*V_inf (O20 strict bound);
beta0 = V_/c_  ->  sharp beta0 = (1-x^2)/(9x^3), x=cos(th/2).

Run:  py -3.12 -u w7/p15_window.py   (deterministic quadrature, no sampling)
"""
import numpy as np
from o21_curv_linear import moments, slice_quad

MODE = "3x"  # label only


def g(x):
    """beta0 (sharp, uncorrected) as a function of x=cos(th/2)."""
    return (1.0 - x * x) / (9.0 * x ** 3)


def solve_x(target, lo=1e-3, hi=1.0 - 1e-9):
    """solve g(x)=target on (lo,hi); g is strictly decreasing."""
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(mid) > target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def th_over_pi(x):
    return 2.0 * np.arccos(x) / np.pi


L = []


def P(s=""):
    print(s)
    L.append(s)


THETAS = [0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4]
B25 = 2.5

P("=" * 100)
P("P15  window recalculation:  3x E[H] in O21's c  ->  beta0  ->  theta-window on beta in [2.5,6]")
P("=" * 100)

# ---------------------------------------------------------------- [C1] reproduce O21's c(2.5)
P("\n[C1] c(2.5) recomputed from O21's model (independent run of the same quadrature)")
P("   theta      E[Q|0]=E[H]/beta     E[Phi''|0]      c(2.5)          O21 md table    rel.diff")
REF = {0.0: 0.291970, np.pi / 4: 0.273852, np.pi / 2: 0.223273, 3 * np.pi / 4: 0.152250}
rows = {}
for th in THETAS:
    m = moments(B25, th)
    c = m["EB"] / B25
    rows[th] = (m["Eq"], m["EPhi_t"], c)
    P(f"   {th:7.4f}   {m['Eq']: .10e}   {m['EPhi_t']: .10f}   {c: .10e}   {REF[th]:.6f}"
      f"     {abs(c-REF[th])/REF[th]:.1e}")

# ---------------------------------------------------------------- [C2] three c variants
P("\n[C2] the two c-sets (+ the exact decomposition of a 3x error in E[H])")
P("   variant A: uncorrected          c_A = E[Q] + E[Phi'']/beta           (O21 as published)")
P("   variant B: task-literal  /3     c_B = c_A/3                          (whole c divided by 3)")
P("   variant C: exact propagation    c_C = E[Q]/3 + E[Phi'']/beta         (only H_vv carries the 3x)")
P("   theta      c_A            c_B=c_A/3      c_C (exact)     c_A/c_C")
cB, cC = {}, {}
for th in THETAS:
    Eq, Ep, cA = rows[th]
    cB[th] = cA / 3.0
    cC[th] = Eq / 3.0 + Ep / B25
    P(f"   {th:7.4f}   {cA: .10f}   {cB[th]: .10f}   {cC[th]: .10f}   {cA/cC[th]: .4f}")

# ---------------------------------------------------------------- [C3] c_inf, and the 3x at c_inf is exact
P("\n[C3] c_inf = cos(th/2)/(2N): uncorrected vs 3x-corrected (Phi''/beta -> 0 at beta->inf,")
P("     so the 3x on c_inf is EXACT: c_inf -> c_inf/3, i.e. beta0 -> 3*beta0)")
P("   theta     c_inf=cos/4     c_inf/3        c(200) meas.    c(200)/c_inf   c_corr(200)   (c_inf/3)")
for th in THETAS:
    ci = 0.25 * np.cos(th / 2.0)
    m = moments(200.0, th)
    c200 = m["EB"] / 200.0
    c200c = m["Eq"] / 3.0 + m["EPhi_t"] / 200.0
    P(f"   {th:7.4f}   {ci: .10f}   {ci/3: .10f}   {c200: .10f}   {c200/ci: .6f}"
      f"     {c200c: .10f}   {ci/3: .10f}")

# ---------------------------------------------------------------- [C4] beta0, four tables
P("\n[C4] beta0 = V_inf/c_inf,  V_inf = tan^2(th/2)/36")
P("     sharp : Var <= V_inf          -> beta0 = (1-x^2)/(9x^3)")
P("     strict: Var <= 6.51 V_inf     -> beta0 = 6.51*(1-x^2)/(9x^3)")
P("     corrected (3x on c_inf)       -> multiply either by 3")
P("   theta     V_inf        c_inf       b0_sharp   b0_sharp_x3   b0_strict  b0_strict_x3")
for th in THETAS:
    x = np.cos(th / 2.0)
    Vi = np.tan(th / 2.0) ** 2 / 36.0
    ci = 0.25 * x
    b0 = Vi / ci
    P(f"   {th:7.4f}   {Vi: .7f}   {ci: .7f}   {b0:9.6f}   {3*b0:11.6f}   {6.51*b0:9.6f}   {3*6.51*b0:11.6f}")
P("   cross-check vs O22 md table (sharp): 0.02063 / 0.15713 / 1.69227  -> reproduced above")

# ---------------------------------------------------------------- [C5] windows
P("\n[C5] theta-windows on beta in [2.5,6]  (theta in units of pi; x=cos(th/2))")
P("     full = {beta0<=2.5}: whole [2.5,6] closed;  partial = {2.5<beta0<=6}: only [beta0,6];")
P("     none = {beta0>6}")
P("   bound              K      x(b0=2.5)  th/pi(b0=2.5)   x(b0=6)   th/pi(b0=6)   verdict")
for name, K in [("sharp", 1.0), ("strict", 6.51), ("sharp x3", 3.0), ("strict x3", 3.0 * 6.51)]:
    x25 = solve_x(2.5 / K)
    x6 = solve_x(6.0 / K)
    P(f"   {name:16s} {K:8.2f}   {x25:.6f}     {th_over_pi(x25):.5f}       {x6:.6f}    "
      f"{th_over_pi(x6):.5f}      full: th<= {th_over_pi(x25):.4f}pi ; "
      f"partial: ({th_over_pi(x25):.4f},{th_over_pi(x6):.4f}]pi ; none: > {th_over_pi(x6):.4f}pi")

# ---------------------------------------------------------------- [C6] provenance of the 3x
P("\n[C6] provenance check: is O21's c on the with-1/N branch or the inflated (no-1/N) branch?")
P("     O15 log [V3]/[V6] and H1A archive (all at N=2, theta=0, beta=2.5):")
P("        with 1/N  : H1A archived 0.553061 ; O15 MC 0.553129 ; O21 model 0.553129 (see [C1] EQ*beta)")
P("        no 1/N    : 2 x that = 1.1062  (would be the inflated branch; contradicted by 3 independent numbers)")
P("     O15 log [V6] (N=3, Z_N centre, beta=6): with 1/N = 0.1021 ; BETASTAR-N Sec.4 = 0.3060")
P(f"        ratio 0.3060/0.1021 = {0.3060/0.1021:.4f}  ~ N=3  -> the 3x lives in BETASTAR-N's N=3 column,")
P("        NOT in O21's N=2 c (O21 has no free N factor: S=beta(1-ReTr/N) with N=2 already).")
P("     general-N coefficient (exact, rho->0): c_inf = cos(th/2)/(2N) = Re tr U_p/(4N)  [T]")
P("        N=2 -> cos/4 (O21) ;  N=3 -> cos/6 (1.5x smaller, i.e. beta0 x1.5 at fixed V_inf)")

# ---------------------------------------------------------------- [C7] consistency of the corrected window
P("\n[C7] does the 3x-corrected window still contain theta<0.5pi?  (decision line (a)/(b))")
for name, K in [("sharp x3", 3.0), ("strict x3", 3.0 * 6.51)]:
    x25 = solve_x(2.5 / K)
    P(f"   {name:10s}: th_c(beta0=2.5) = {th_over_pi(x25):.5f} pi  -> "
      f"{'NON-EMPTY, contains 0.5pi' if th_over_pi(x25) >= 0.5 else 'EMPTY at 0.5pi'}")
P("   (b) would need th_c<0 for every theta, i.e. beta0>2.5 even at theta->0: false, beta0(theta=0)=0 identically.")

# ---------------------------------------------------------------- [C8] quadrature convergence
P("\n[C8] convergence (npan 240 -> 480, beta=2.5)")
import o21_curv_linear as o21
keep = (o21.RHO, o21.DRHO)
for th in [0.0, np.pi / 2]:
    m1 = o21.moments(2.5, th)
    o21.RHO, o21.DRHO = slice_quad(npan=480, nnode=40)
    m2 = o21.moments(2.5, th)
    P(f"   theta={th:6.4f}  c: {m1['EB']/2.5:.12f} -> {m2['EB']/2.5:.12f}")
    o21.RHO, o21.DRHO = keep

with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as f:
    f.write("\n".join(L) + "\n")
print("\n[log written]")
