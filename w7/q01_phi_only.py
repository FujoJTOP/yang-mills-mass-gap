"""
Q01 - H1 on the Phi-only measure (S_W = 0, beta = 0): compute G''(0), Var(Phi'), E(Phi'')
       and the S-D split.   py -3.12, numpy only, deterministic quadrature (no sampling).

Conventions, verbatim from the corpus (not from memory):
  flat chart      dU = J dxi ,  J = sin^2(rho/2)/rho^2 ,  Phi = -log J        (sec03-failure.tex eq:model1)
  Phi             = sum_{i<j} g(lam_i-lam_j), g(d) = -2 log(sin(d/2)/(d/2))
                  = h(rho),  rho = |xi| ,  h(rho) = -2 log(sin(rho/2)/(rho/2))  (w7/p16_phiquant.py)
  measure         dmu = e^{-S} J d^3xi ; beta = 0  =>  dmu = J d^3xi = Haar(SU(2))
  H1              Var_mu(<v,grad S>|0) <= E_mu[H_vv|0] + W_full - delta,  V = S_W + Phi
  G_v = <v,grad S_W>,  H_vv = <v,Hess S_W v>,  Phi' = <v,grad Phi>,  Phi'' = <v,Hess Phi v>
  W_geo = E[Phi''] - Var(Phi')   (ch05-theorems.tex eq:wgeo; sec02-scalar.tex table)
  conditioning:  v = e_3 (unit),  a_0 = <v,xi> = xi_3 ,  slice xi_3 = 0

On the Phi-only measure: S_W = 0 => G_v == 0, H_vv == 0, V = Phi, total force = Phi'.
Quadrature: rho = 2pi(1-w^2) clusters nodes at the cut locus rho = 2pi (antipode -I), where
Phi' ~ 2/(2pi-rho), Phi'' ~ 2/(2pi-rho)^2, J ~ (2pi-rho)^2/4 : every product is smooth in w.
"""
import numpy as np

LOG = open(r"D:\Dev\math\w7\q01_phi_only.log", "w", encoding="utf-8")


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.write(s + "\n")
    LOG.flush()


TWO_PI = 2.0 * np.pi
N = 2          # SU(2)
FLOOR = N / 12.0


def _small(r):
    return r < 0.05


def hp(r):
    """h'(rho) = 2/rho - cot(rho/2)."""
    r = np.asarray(r, float)
    sm = _small(r)
    out = np.empty_like(r)
    out[sm] = r[sm] / 6 + r[sm] ** 3 / 360 + r[sm] ** 5 / 15120
    out[~sm] = 2 / r[~sm] - 1 / np.tan(r[~sm] / 2)
    return out


def hpp(r):
    """h''(rho) = csc^2(rho/2)/2 - 2/rho^2."""
    r = np.asarray(r, float)
    sm = _small(r)
    out = np.empty_like(r)
    out[sm] = 1 / 6 + r[sm] ** 2 / 120 + r[sm] ** 4 / 3024
    out[~sm] = 0.5 / np.sin(r[~sm] / 2) ** 2 - 2 / r[~sm] ** 2
    return out


def div(a, b, lim):
    """a/b with the r -> 0 limit `lim` filled in (no 0/0)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    return np.divide(a, b, out=np.full(np.broadcast(a, b).shape, float(lim)), where=(b > 1e-12))


def rat(r):
    """h'(rho)/rho  -> 1/6 at 0."""
    r = np.asarray(r, float)
    return div(hp(np.maximum(r, 1e-9)), np.maximum(r, 1e-9), 1.0 / 6)


def simp(g, a=0.0, b=1.0, n=200000):
    """composite Simpson on [a,b] with n (even) panels."""
    x = np.linspace(a, b, n + 1)
    y = g(x)
    h = (b - a) / n
    return float(h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum()))


def quad(f, lo_w=0.0, n=200000):
    """int_0^{2pi} f(rho) drho  (optionally cut at rho = 2pi - eps via lo_w = sqrt(eps/2pi))."""
    lo = max(lo_w, 0.0)

    def g(w):
        rho = TWO_PI * (1 - w ** 2)
        return f(rho) * (4 * np.pi * w)           # |drho/dw| = 4 pi w

    return simp(g, lo, 1.0, n)


def Cin(x, n=20000):
    def g(t):
        return np.divide(1 - np.cos(t), t, out=np.zeros_like(t), where=(t > 1e-12))
    return simp(g, 0.0, float(x), n)


say("=" * 100)
say("Q01  Phi-only (S_W = 0, beta = 0, SU(2), single plaquette, flat chart, v = e_3, a_0 = 0)")
say("=" * 100)
Z3 = 4 * np.pi * quad(lambda r: np.sin(r / 2) ** 2)
say(f"[0] normalisation:  4pi*int sin^2(rho/2)drho = {Z3:.10f}   (4pi^2 = {4*np.pi**2:.10f})  <- quadrature check")

say("=" * 100)
say("[1] UNCONDITIONED (marginal over a): Stein/IBP  E[Phi''] = E[(Phi')^2]  =>  Var(Phi') = E(Phi'') - (E Phi')^2")
E_pp = quad(lambda r: np.sin(r / 2) ** 2 * ((1 / 3) * hpp(r) + (2 / 3) * rat(r))) / np.pi
E_p2 = quad(lambda r: np.sin(r / 2) ** 2 * ((1 / 3) * hp(r) ** 2)) / np.pi
say(f"    E[Phi''_vv]        = {E_pp:.12f}")
say(f"    E[(Phi'_v)^2]      = {E_p2:.12f}")
say(f"    |diff|/E          = {abs(E_pp - E_p2)/E_pp:.3e}     (IBP; E[Phi'_v] = 0 by xi -> -xi)")
say(f"    => Var_uncond(Phi'_v) = E(Phi'') = {E_p2:.12f}  EXACTLY: slack = (E Phi')^2 = 0")

say("=" * 100)
say("[2] CONDITIONED at a_0 = 0 (the object H1 uses): slice xi_3 = 0")
say("    dmu|slice = J dxi_1 dxi_2 = sin^2(r/2)/r dr dphi ; Phi radial => d_3 Phi = h'(rho) xi_3/rho = 0 on slice")
Fs = 0.5 * Cin(TWO_PI)
say(f"    int sin^2(r/2)/r dr = (1/2)Cin(2pi) = {Fs:.12f}    Cin(2pi) = {2*Fs:.12f}")
E_pp0 = quad(lambda r: rat(r) * div(np.sin(r / 2) ** 2, r, 0.0)) / Fs
say(f"    Var(Phi'|0)        = {0.0:.12f}    (exact: Phi' == 0 identically on the slice)")
say(f"    Phi''|0            = h'(r)/r >= 1/6 pointwise")
say(f"    E[Phi''|0]         = {E_pp0:.12f}")
say(f"    1/(2 Cin(2pi))     = {1/(4*Fs):.12f}    <- corpus W_0 = 0.20511530 (H1A-EXACT S2.1; p16 log [6])")
say(f"    N/12               = {FLOOR:.12f}    ratio = {E_pp0/FLOOR:.6f}")

say("=" * 100)
say("[3] EXACT marginal curvature: p(a) = (2pi/Z) int_{|a|}^{2pi} sin^2(rho/2)/rho drho")
say("    G''_Phi(a) = -d^2/da^2 log p = f'(a)/(Fs-F(a)) - [f(a)/(Fs-F(a))]^2 ,  f = sin^2(a/2)/a, F' = f")
say("    (== E[Phi''|a] - Var(Phi'|a) by the conditioning identity; must reproduce [2] at a = 0)")


def f_of(a):
    a = np.atleast_1d(np.asarray(a, float))
    out = np.empty_like(a)
    sm = a < 1e-3
    out[sm] = a[sm] / 4 - a[sm] ** 3 / 48
    out[~sm] = np.sin(a[~sm] / 2) ** 2 / a[~sm]
    return out


def fp_of(a):
    a = np.atleast_1d(np.asarray(a, float))
    out = np.empty_like(a)
    sm = a < 1e-3
    out[sm] = 0.25 - a[sm] ** 2 / 16 + a[sm] ** 4 / 480
    out[~sm] = np.sin(a[~sm] / 2) * (a[~sm] * np.cos(a[~sm] / 2) - np.sin(a[~sm] / 2)) / a[~sm] ** 2
    return out


ag = np.concatenate([[0.0], np.geomspace(1e-3, TWO_PI - 1e-3, 160)])
G2 = np.empty(ag.size)
for i, a in enumerate(ag):
    rem = Fs - 0.5 * Cin(float(a), 4000)
    G2[i] = (fp_of(a)[0] / rem) - (f_of(a)[0] / rem) ** 2
say(f"    G''_Phi(0)         = {G2[0]:.12f}    (analytic f'(0)/Fs = {0.25/Fs:.12f} = 1/(2Cin(2pi)) )")
say(f"    min over a-grid    = {G2.min():.12f}    at a = {ag[G2.argmin()]:.6f}")
say("    a = 0.5,1,2,3,4,5,6: " + " ".join(f"{G2[np.argmin(abs(ag-x))]:.5f}" for x in (0.5, 1, 2, 3, 4, 5, 6)))
pos = np.where(G2 > 0)[0]
across = ag[pos[-1]] if pos.size < ag.size else float("nan")
say(f"    zero crossing      : G''_Phi > 0 up to a = {across:.4f}, negative beyond")
say(f"    edge               : G''_Phi(2pi-1e-3) = {G2[-1]:.6e}   (p ~ (2pi-a)^3 => -d^2 log p ~ -3/(2pi-a)^2 -> -infinity)")
say(f"    => UNIFORM in a FAILS on the flat chart: G''_Phi(0) = {G2[0]:.6f} >= N/12 but inf_a G''_Phi = -inf at the cut locus")
say(f"    S-D at a=0: S = E[Phi''|0] = {E_pp0:.6f}, D = Var(Phi'|0) = 0")

say("=" * 100)
say("[4] BL leg: Hess Phi >= (N/12)I => Var(g) <= (12/N)E|grad g|^2 ; g = Phi'_v, grad g = Hess Phi . v")
say("    |Hess Phi e_3|^2 = h''^2 cos^2 + (h'/rho)^2 sin^2 ; angular average (1/3, 2/3)")
say("    cut at rho = 2pi - eps and let eps -> 0:")
for eps in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5):
    lo = np.sqrt(eps / TWO_PI)
    num = quad(lambda r: np.sin(r / 2) ** 2 * ((1 / 3) * hpp(r) ** 2 + (2 / 3) * rat(r) ** 2), lo_w=lo) / np.pi
    say(f"    eps={eps:8.0e}:  E_eps|grad Phi'|^2 = {num:12.6e}   BL bound (12/N)E = {6*num:12.6e}   eps*B = {eps*6*num:.6f}")
say(f"    measured Var(Phi'_v) = {E_p2:.6f}  ->  the BL bound with f = Phi' DIVERGES (vacuous), it is not loose.")

say("=" * 100)
say("[5] LEDGER (v = e_3 unit; L links iid => sum w_l^2 = 1 gives the same number for every L, incl. L=2)")
say("    G_v = d_v S_W        = 0          [exact: S_W == 0]")
say("    H_vv = v^T Hess S_W v = 0         [exact]")
say("    total force = Phi' ;  Phi'' = h''cos^2 + (h'/rho)sin^2 ;  min Phi'' = N/12 = %.6f" % FLOOR)
say(f"    conditioned  : Var(Phi'|0) = 0 <= E[Phi''|0] = {E_pp0:.7f}   margin = {E_pp0:.7f}")
say(f"    unconditioned: Var(Phi')  = E(Phi'') = {E_p2:.7f}          margin = 0 exactly")
say(f"    G''_Phi(0) = E[Phi''|0] - Var(Phi'|0) = {E_pp0:.7f}")
LOG.close()
