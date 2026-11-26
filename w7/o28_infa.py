"""
O28 - independent recomputation of inf_a beta*(a) for the single-plaquette SU(2) toy.

Model (verbatim from research/v3/species/H1A-EXACT.md 2.1, = O09):
  flat chart xi in R^3 = su(2), rho=|xi|, u=xi_3/rho, J = sin^2(rho/2)/rho^2, Phi = -log J
  S = beta (1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) u),  v = e_3, slice xi_3 = a
  V = S + Phi ;  G(a) = -log[ int_slice e^{-V} s ds ]  (s = transverse radius)
  G''(a) = E[V_33 | a] - Var(V_3 | a) = E[H_vv] + E[Phi''] - Var(force)

INDEPENDENCE vs O09: O09 gets V_3, V_33 by 5-point finite differences in x_3 (h=1e-4).
Here both are the EXACT analytic coordinate derivatives (derived by chain rule below),
so the two routes share only the model, not the differentiation.

Second independent route (method A): G(a) itself by quadrature, then 5-point stencil
in a on G -> G''. Compared against the analytic route (method B) at every reported a.
"""
import sys
import numpy as np

LOG = open(r"D:\Dev\math\w7\o28_infa.log", "w", encoding="utf-8")


def say(*x):
    line = " ".join(str(t) for t in x)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


np.set_printoptions(precision=10, suppress=True)
TWOPI = 2.0 * np.pi
GL = {}


def gl(n):
    if n not in GL:
        x, w = np.polynomial.legendre.leggauss(n)
        GL[n] = (0.5 * (x + 1.0), 0.5 * w)
    return GL[n]


def nodes(a, n):
    """s in [0, s_max] on the slice xi_3=a, s_max = sqrt((2pi)^2 - a^2).  Exact domain
    (rho < 2pi); the weight carries J -> 0 at rho=2pi so the integrand vanishes there."""
    smax = np.sqrt(max(TWOPI ** 2 - a * a, 0.0))
    nd, wt = gl(n)
    return smax * nd, smax * wt


def pieces(rho, u, a, beta, th):
    """S, Phi and their x_3-derivatives.  A(rho)=-beta c cos(rho/2), B(rho)=beta s sin(rho/2)u."""
    c, s2 = np.cos(th / 2.0), np.sin(th / 2.0)
    h = rho / 2.0
    sh, ch = np.sin(h), np.cos(h)
    cth = np.cos(h)
    sth = sh
    A_p = 0.5 * beta * c * sth
    A_pp = 0.25 * beta * c * cth
    Bp = 0.5 * beta * s2 * cth
    Bpp = -0.25 * beta * s2 * sth
    B = beta * s2 * sth
    S = beta * (1.0 - cth * c + sth * s2 * u)
    Phi = -2.0 * np.log(np.maximum(sh / np.maximum(h, 1e-300), 1e-300))
    Ph_p = 2.0 / rho - cth / np.maximum(sh, 1e-300)
    Ph_pp = -2.0 / rho ** 2 + 0.5 / np.maximum(sh, 1e-300) ** 2
    # d/dxi_3 and d^2/dxi_3^2 at fixed transverse radius s
    S_3 = (A_p + Bp * u) * u + B * (1.0 - u * u) / rho
    Phi_3 = Ph_p * u
    S_33 = ((A_pp + Bpp * u) * u * u
            + 2.0 * Bp * u * (1.0 - u * u) / rho
            + (A_p + Bp * u) * (1.0 - u * u) / rho
            - 3.0 * B * u * (1.0 - u * u) / rho ** 2)
    Phi_33 = Ph_pp * u * u + Ph_p * (1.0 - u * u) / rho
    return S, Phi, S_3, Phi_3, S_33, Phi_33


def slice(a, beta, th, n=1200):
    s, wt = nodes(a, n)
    rho = np.sqrt(s * s + a * a)
    rho = np.maximum(rho, 1e-14)
    u = a / rho
    S, Phi, S3, P3, S33, P33 = pieces(rho, u, a, beta, th)
    V = S + Phi
    w = np.exp(-(V - V.min())) * s * wt
    Z = w.sum()
    w = w / Z
    V3 = S3 + P3
    V33 = S33 + P33
    m3 = float((w * V3).sum())
    var = float((w * (V3 - m3) ** 2).sum())
    out = dict(Z=Z, V=V,
               Phi2=float((w * P33).sum()), Hvv=float((w * S33).sum()),
               varF=var, Gpp=float((w * V33).sum()) - var)
    return out


def Gof(a, beta, th, n=1200):
    """-log Z(a) with the (dominant) e^{-V_min} factor kept: use log-sum-exp."""
    s, wt = nodes(a, n)
    rho = np.sqrt(s * s + a * a)
    rho = np.maximum(rho, 1e-14)
    u = a / rho
    S, Phi, *_ = pieces(rho, u, a, beta, th)
    V = S + Phi
    vmin = V.min()
    Z = float((np.exp(-(V - vmin)) * s * wt).sum())
    return vmin - np.log(Z)


def GppA(a, beta, th, h=2e-3, n=1200):
    g = lambda z: Gof(z, beta, th, n)
    return (-g(a + 2 * h) + 16 * g(a + h) - 30 * g(a) + 16 * g(a - h) - g(a - 2 * h)) / (12 * h * h)


say("=" * 96)
say("O28  independent recomputation of inf_a beta*(a);  model = H1A-EXACT 2.1 (single plaquette SU(2))")
say(f"numpy {np.__version__}  python {sys.version.split()[0]}  GL nodes n=1200 (convergence check vs 600/2400 below)")
say("=" * 96)
say("")
say("[0] VALIDATION at a=0 against H1A-EXACT 2.1 / O09 [4]")
vals = ((np.pi, "theta=pi", (0.2051153, 0.0, 0.069121, 0.135995)),
        (0.0, "theta=0 ", (0.1768636, 0.553061, 0.0, 0.7299241)))
for th, tag, ref in vals:
    r = slice(0.0, 2.5, th)
    say(f"  {tag}: E[Phi'']={r['Phi2']:.7f} (ref {ref[0]:.7f})  E[H_vv]={r['Hvv']:+.3e} (ref {ref[1]:.3f})"
        f"  Var(force)={r['varF']:.7f} (ref {ref[2]:.7f})  G''={r['Gpp']:+.7f} (ref {ref[3]:.7f})")
say("")

say("[0b] quadrature convergence at a in {0, 0.709, 1.9} (beta=2.5, theta=pi): G'' for n = 600/1200/2400")
for a in (0.0, 0.70906, 1.9):
    row = [slice(a, 2.5, np.pi, n)["Gpp"] for n in (600, 1200, 2400)]
    say(f"  a={a:6.5f}  " + "  ".join(f"{v:+.12f}" for v in row) + f"   drift={row[2]-row[1]:+.1e}")
say("")

say("[1] two independent routes to G''(a) (beta=2.5, theta=pi): analytic vs G(a)+stencil")
say(f"  {'a':>8} {'Gpp analytic(B)':>18} {'Gpp stencil(A)':>18} {'A-B':>11}")
for a in (0.0, 0.25, 0.5, 0.70906, 0.75, 1.0, 1.5, 1.9, 2.0, 2.5, 2.97192, 3.0):
    b = slice(a, 2.5, np.pi)["Gpp"]
    am = GppA(a, 2.5, np.pi)
    say(f"  {a:8.5f} {b:+18.10f} {am:+18.10f} {am-b:+11.2e}")
say("")

say("[2] G''(a) curve at beta=2.5, theta=pi : 401 points on a in [0, pi]")
AG = np.linspace(0.0, np.pi, 401)
GP = np.array([slice(float(a), 2.5, np.pi)["Gpp"] for a in AG])
say(f"  grid: {len(AG)} points, da = {AG[1]-AG[0]:.6f}")
say(f"  G''(0) = {GP[0]:+.10f}   G''(pi) = {GP[-1]:+.10f}   min = {GP.min():+.10f} at a = {AG[GP.argmin()]:.5f}"
    f"   max = {GP.max():+.10f} at a = {AG[GP.argmax()]:.5f}")
neg = AG[GP < 0]
say(f"  sign changes: {int(np.sum(np.diff(np.sign(GP)) != 0))};  #grid points with G''<0 = {len(neg)}"
    f"  (a-range of negativity: {('%.5f..%.5f' % (neg[0], neg[-1])) if len(neg) else 'none'})")
say("  curve (every 10th point):")
for i in range(0, len(AG), 10):
    say(f"    a={AG[i]:.4f}  G''={GP[i]:+.8f}")
say("")

say("[3] zeros of G''(a) at beta=2.5 by bisection on the 401-pt sign changes")
def bis_a(lo, hi, beta=2.5, th=np.pi, it=80):
    f = lambda z: slice(z, beta, th)["Gpp"]
    flo = f(lo)
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return 0.5 * (lo + hi)
roots = []
for i in range(len(AG) - 1):
    if GP[i] * GP[i + 1] < 0:
        roots.append(bis_a(AG[i], AG[i + 1]))
for r in roots:
    say(f"  root a* = {r:.6f}   G''(a*-1e-6)={slice(r-1e-6,2.5,np.pi)['Gpp']:+.2e}"
        f"  G''(a*+1e-6)={slice(r+1e-6,2.5,np.pi)['Gpp']:+.2e}")
say("")

say("[4] beta*(a): bisection in beta on G''(beta;a)=0 (theta=pi)")
def beta_star(a, lo=1e-6, hi=40.0, it=60):
    f = lambda b: slice(a, b, np.pi)["Gpp"]
    flo, fhi = f(lo), f(hi)
    if not (flo > 0 > fhi):
        return None, flo, fhi
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), flo, fhi

say(f"  {'a':>7} {'beta*(a)':>12} {'G''(beta=2.5)':>15}")
tab = {}
for a in (0.0, 0.5, 1.0, 1.5, 1.9, 2.0, 2.5):
    b, flo, fhi = beta_star(a)
    tab[a] = b
    say(f"  {a:7.3f} {('%.5f' % b) if b else 'none':>12} {slice(a,2.5,np.pi)['Gpp']:+15.8f}")
say("")

say("[5] inf over a: beta*(a) on a fine grid, unconstrained minimisation around the minimum")
AF = np.linspace(0.0, np.pi, 161)
BS = np.array([beta_star(float(a))[0] for a in AF])
i = int(np.nanargmin(BS))
say(f"  coarse scan {len(AF)} pts: min beta* = {BS[i]:.5f} at a = {AF[i]:.3f}")
lo, hi = AF[max(i - 1, 0)], AF[min(i + 1, len(AF) - 1)]
gr = np.linspace(lo, hi, 61)
bg = np.array([beta_star(float(a))[0] for a in gr])
j = int(np.nanargmin(bg))
say(f"  refined scan {len(gr)} pts on [{lo:.3f},{hi:.3f}]: min beta* = {bg[j]:.6f} at a = {gr[j]:.5f}")
# golden-section on beta*(a)
phi = (np.sqrt(5) - 1) / 2
x1, x2 = gr[max(j - 1, 0)], gr[min(j + 1, len(gr) - 1)]
for _ in range(60):
    m1 = x2 - phi * (x2 - x1)
    m2 = x1 + phi * (x2 - x1)
    if beta_star(float(m1))[0] < beta_star(float(m2))[0]:
        x2 = m2
    else:
        x1 = m1
amin = 0.5 * (x1 + x2)
bmin = beta_star(amin)[0]
say(f"  golden-section: a_min = {amin:.6f}   beta*(a_min) = {bmin:.6f}")
say(f"  => inf_a beta*(a) = {bmin:.5f}")
say("")
say("[6] judgement inputs")
say(f"  beta*(0)             = {tab[0.0]:.5f}")
say(f"  inf_a beta*(a)       = {bmin:.5f}")
say(f"  beta_c(SU2)          = 2.2986")
say(f"  inf < beta_c ?       = {bmin < 2.2986}")
say(f"  beta*(0) > beta_c ?  = {tab[0.0] > 2.2986}")
say(f"  G''(a) >= 0 for all a in [0,pi] at beta=2.5 ? = {bool(np.all(GP >= 0))}")
say("")
say("[7b] cross-check of each root: analytic and stencil routes evaluated AT beta*(a)")
for a in (0.0, 0.5, 1.0, 1.5, 1.9, 2.0, 2.5):
    b = tab[a]
    if b is None:
        continue
    say(f"  a={a:5.2f}  beta*={b:.6f}  G''(B)={slice(a,b,np.pi)['Gpp']:+.3e}  G''(A)={GppA(a,b,np.pi):+.3e}")
say("")
say("[7c] is the O09-vs-O28 gap of ~2e-4 in beta* caused by the RMAX=2pi-1e-3 truncation?")
_saved = None
for a in (1.0, 1.9):
    r12 = slice(a, tab[a], np.pi, n=2400)["Gpp"]
    s_, wt_ = nodes(a, 2400)
    # same quadrature but with the domain shaved to 2pi-1e-3 exactly as O09 does
    smax = np.sqrt((TWOPI - 1e-3) ** 2 - a * a)
    nd, wt = gl(2400)
    s2_ = smax * nd
    rho = np.sqrt(s2_ ** 2 + a * a)
    u = a / rho
    S, Phi, S3, P3, S33, P33 = pieces(rho, u, a, tab[a], np.pi)
    V = S + Phi
    w = np.exp(-(V - V.min())) * s2_ * smax * wt
    w = w / w.sum()
    V3 = S3 + P3
    m3 = float((w * V3).sum())
    g_shaved = float((w * (S33 + P33)).sum()) - float((w * (V3 - m3) ** 2).sum())
    say(f"  a={a:5.2f} beta={tab[a]:.6f}: G''(domain 2pi)={r12:+.3e}   G''(domain 2pi-1e-3)={g_shaved:+.3e}"
        f"   shift={g_shaved-r12:+.2e}")
say("")
say("[8] corpus survivor reading (V3-COUNTEREXAMPLE-SCOPE.md:38-47): evaluate at the coarse")
say("    minimiser a0* of S_eff(a) = -log Z(a), not uniformly in a0.  Scan a in (-pi,pi).")
AALL = np.linspace(-np.pi + 1e-6, np.pi - 1e-6, 401)
for beta in (2.5, 2.9, 3.5):
    Seff = np.array([Gof(float(a), beta, np.pi, 800) for a in AALL])
    k = int(np.argmin(Seff))
    astar = float(AALL[k])
    say(f"  beta={beta}: arg min S_eff = a0* = {astar:+.4f}   G''(a0*) = {slice(astar,beta,np.pi)['Gpp']:+.6f}"
        f"   G''(0) = {slice(0.0,beta,np.pi)['Gpp']:+.6f}")
say("")
say("[9] beta* at the coarse minimiser a0*(beta) (self-consistent: a0* depends on beta)")
def astar_of(beta):
    S = np.array([Gof(float(a), beta, np.pi, 800) for a in AALL])
    i = int(np.argmin(S))
    lo, hi = AALL[max(i - 1, 0)], AALL[min(i + 1, len(AALL) - 1)]
    gr = np.linspace(lo, hi, 41)
    Sg = np.array([Gof(float(a), beta, np.pi, 800) for a in gr])
    return float(gr[int(np.argmin(Sg))])

def Gpp_at_min(beta, it=6):
    a = astar_of(beta)
    return a, slice(a, beta, np.pi)["Gpp"]

for beta in (2.0, 2.5, 3.0, 4.0, 5.0):
    a, g = Gpp_at_min(beta)
    say(f"  beta={beta:.1f}: a0*={a:+.4f}  G''(a0*)={g:+.6f}")
lo, hi = 1.0, 6.0
for _ in range(40):
    mid = 0.5 * (lo + hi)
    if Gpp_at_min(mid)[1] > 0:
        lo = mid
    else:
        hi = mid
say(f"  beta* at the coarse minimiser = {0.5*(lo+hi):.5f}   (a0* there = {Gpp_at_min(0.5*(lo+hi))[0]:+.4f})")
say("")
say("[7] where the negativity lives (beta=2.5): a-window and gap to the vacuum")
say("  vacuum of this slice family sits at xi_3 = -pi (theta=pi).  conditioning point xi_3=a.")
for a in (0.0, 0.70906, 1.0, 1.9, 2.97192):
    say(f"    a={a:7.5f}  |a - (-pi)| = {abs(a + np.pi):.4f}   G''={slice(a,2.5,np.pi)['Gpp']:+.6f}")
say("=" * 96)
LOG.close()
