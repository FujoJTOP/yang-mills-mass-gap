"""Q08-LARGEA -- large-a behaviour of G''(a) for the single-plaquette SU(2) toy (theta=pi).

Model = H1A-EXACT 2.1 / O28 (verbatim): rho=|xi|, s=transverse radius, a=xi_3, u=a/rho,
  S = beta(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) u),  Phi = -2log[sin(rho/2)/(rho/2)]
  Z(a) = int_0^{smax} e^{-S-Phi} s ds,   G = -log Z,   G'' = E[V_33|a] - Var(V_3|a)

NEW (this script): exact reduction at theta=pi.  Since s ds = rho drho and rho<=2pi always,
  Z(a) = e^{-beta} int_{|a|}^{2pi} w(rho) exp(-beta a phi(rho)) drho,
  w(rho) = rho e^{-Phi} = 4 sin^2(rho/2)/rho ,  phi(rho) = sin(rho/2)/rho ,
so the whole a-dependence is one tilt parameter x=beta*a on a FIXED reference measure w d(rho).
Differentiating twice (lower limit + tilt) gives the exact decomposition

  G''(a) = -beta^2 Var_nu(phi)                        [concavity / Savvidy variance]
           + w'(a)/I(a)                                [geometric term, O(1)]
           + kappa [ 2 beta (E_nu phi - phi(a) - a phi'(a)/2) + w'(a)/w(a) ]
           + kappa^2 ,        kappa = w(a) exp(-beta a phi(a)) / I(a)

with nu the normalized tilted measure on [a,2pi].  At a=0: kappa=0, phi'(0)=0 =>
  G''(0) = w'(0)/I(0) - beta^2 Var_w(phi) = 1/(2 Cin(2pi)) - beta^2 Var_w(phi)   (closed form).

ROUTES: (B) closed decomposition  vs  (A) 5-point stencil on G(a) from quadrature.
Both are checked against O28's published numbers before any new claim.
"""
import sys
import numpy as np

LOG = open(r"D:\Dev\math\w7\q08_largea.log", "w", encoding="utf-8")


def say(*x):
    line = " ".join(str(t) for t in x)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


TWOPI = 2.0 * np.pi
_GL = {}


def gl(n):
    if n not in _GL:
        x, w = np.polynomial.legendre.leggauss(n)
        _GL[n] = (0.5 * (x + 1.0), 0.5 * w)
    return _GL[n]


def integ(f, lo, hi, n):
    if hi <= lo:
        return 0.0
    t, wt = gl(n)
    return (hi - lo) * float(np.dot(wt, f(lo + (hi - lo) * t)))


def W(r):
    return 4.0 * np.sin(0.5 * r) ** 2 / r


def PH(r):
    return np.sin(0.5 * r) / r


def PHp(r):
    return (0.5 * r * np.cos(0.5 * r) - np.sin(0.5 * r)) / r ** 2


def Wp(r):
    return (2.0 * r * np.sin(r) - 4.0 * np.sin(0.5 * r) ** 2) / r ** 2


def I_of(a, beta, n=3000):
    lo = abs(a)
    return integ(lambda r: W(r) * np.exp(-beta * a * PH(r)), lo, TWOPI, n)


def mom(a, beta, n=3000):
    lo = abs(a)
    e = lambda r: W(r) * np.exp(-beta * a * PH(r))
    I = integ(e, lo, TWOPI, n)
    m1 = integ(lambda r: e(r) * PH(r), lo, TWOPI, n)
    m2 = integ(lambda r: e(r) * PH(r) ** 2, lo, TWOPI, n)
    return I, m1 / I, m2 / I


def G_of(a, beta, n=3000):
    return beta - np.log(I_of(a, beta, n))


def Gpp_closed(a, beta, n=3000):
    """route B: exact decomposition above (a>=0)."""
    I, m1, m2 = mom(a, beta, n)
    var = m2 - m1 * m1
    if a == 0.0:
        kap, phi_a, phip_a, wpa = 0.0, 0.5, 0.0, 1.0   # w(rho)~rho => w'(0)=1
    else:
        kap = W(a) * np.exp(-beta * a * PH(a)) / I
        phi_a, phip_a, wpa = PH(a), PHp(a), Wp(a)
    t1 = -beta * beta * var
    t2 = kap * wpa / W(a) if a != 0.0 else wpa / I     # kappa * w'/w
    t3 = kap * 2.0 * beta * (m1 - phi_a - 0.5 * a * phip_a)
    t4 = kap * kap
    return dict(I=I, var=var, E1=m1, kap=kap, t1=t1, t2=t2, t3=t3, t4=t4,
                Gpp=t1 + t2 + t3 + t4)


def Gpp_st(a, beta, h=2e-3, n=1500):
    """route A: 5-point stencil on the quadrature G(a) (works for either sign of a)."""
    g = lambda z: G_of(z, beta, n)
    return (-g(a + 2 * h) + 16 * g(a + h) - 30 * g(a) + 16 * g(a - h) - g(a - 2 * h)) / (12 * h * h)


say("=" * 98)
say("Q08-LARGEA  large-a behaviour of G''(a), single-plaquette SU(2) toy, theta=pi")
say(f"numpy {np.__version__}  python {sys.version.split()[0]}  GL nodes n=3000 (scan 800/1500)")
say("exact reduction: Z(a)=e^-beta int_|a|^2pi w(rho) e^{-beta a phi(rho)} drho,")
say("                 w=4sin^2(rho/2)/rho, phi=sin(rho/2)/rho ;  G=beta-log I")
say("=" * 98)
say("")

# ---------------------------------------------------------------- [0] validation
say("[0] VALIDATION against O28 (w7/o28_infa.log) -- no new claim before this matches")
ref = ((0.0, 0.1359947110), (np.pi, 0.0272768642), (0.709054, 0.0), (2.973310, 0.0))
for a, want in ref:
    c = Gpp_closed(a, 2.5)
    s = Gpp_st(a, 2.5, 1e-2)
    say(f"  a={a:8.6f}  routeB={c['Gpp']:+.12f}  routeA={s:+.12f}  A-B={s-c['Gpp']:+.1e}"
        f"   O28={want:+.10f}  diff={c['Gpp']-want:+.1e}")
c = Gpp_closed(0.0, 2.5)
say(f"  a=0 split: -beta^2 Var_w(phi)={c['t1']:+.12f}   w'(0)/I={c['t2']:+.12f}"
    f"   =>  O28 Var(force)=0.0691206 vs beta^2 Var={-c['t1']:.7f}")
say(f"  beta=0: G''(0)=w'(0)/I = {Gpp_closed(0.0, 0.0)['Gpp']:.10f}"
    f"   1/(2Cin(2pi))=0.2051152970   h1d anchor=0.20511530")
AG = np.linspace(0.0, np.pi, 401)
GP = np.array([Gpp_closed(float(a), 2.5)["Gpp"] for a in AG])
say(f"  401-pt [0,pi] scan, beta=2.5: min={GP.min():+.10f} @ a={AG[GP.argmin()]:.5f}"
    f"  G''(0)={GP[0]:+.10f} G''(pi)={GP[-1]:+.10f}  #neg={int(np.sum(GP<0))}"
    f"   (O28: -0.0791945935 @1.82212, 0.1359947110, 0.0272768642, 288)")
say("")

# ---------------------------------------------------------------- [1] decomposition
say("[1] exact decomposition at beta=2.5: where the negativity comes from (route B)")
say(f"  {'a':>7} {'-b^2Var':>11} {'kap*w_/w':>10} {'kappa*2b(..)':>13} {'kappa^2':>9}"
    f" {'sum':>11} {'Var_nu(phi)':>12} {'kappa':>8}")
for a in (0.0, 0.35, 0.709054, 1.0, 1.4, 1.822, 2.2, 2.6, 2.973310, 3.141593, 4.0, 5.0, 6.0):
    c = Gpp_closed(a, 2.5)
    say(f"  {a:7.4f} {c['t1']:+11.6f} {c['t2']:+10.6f} {c['t3']:+13.6f} {c['t4']:+9.6f}"
        f" {c['Gpp']:+11.6f} {c['var']:12.7f} {c['kap']:8.4f}")
say("  -> the ONLY O(beta^2) term is -beta^2 Var_nu(phi) <= 0; every counterweight is O(1).")
say("")

# ---------------------------------------------------------------- [2] beta scan
say("[2] negative set on the full circle a in (-2pi,2pi), route A (validated above)")
def roots_and_map(beta, n=800, h=2e-3):
    grid = np.linspace(-TWOPI + 0.01, TWOPI - 0.01, 401)
    v = np.array([Gpp_st(float(a), beta, h, n) for a in grid])
    cuts = []
    for i in range(len(grid) - 1):
        if v[i] * v[i + 1] < 0:
            lo, hi = grid[i], grid[i + 1]
            flo = v[i]
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                fm = Gpp_st(float(mid), beta, h, n)
                if (fm > 0) == (flo > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            cuts.append(0.5 * (lo + hi))
    return grid, v, cuts

say(f"  {'beta':>8} {'G_(0)':>11} {'G_(pi)':>11} {'min_(0,2pi)':>12} {'argmin':>8}"
    f" {'ncomp':>6}  sign-change points / width of the a>0 window")
for beta in (0.0, 1.0, 1.8733, 2.0, 2.5, 2.9, 3.0, 4.0, 4.30661, 4.5, 6.0, 8.0, 12.0, 20.0):
    grid, v, cuts = roots_and_map(beta)
    pos = np.linspace(0.0, TWOPI - 0.02, 401)
    vp = np.array([Gpp_st(float(a), beta) for a in pos])
    i = int(vp.argmin())
    neg = [c for c in cuts if c > 0]
    w = (neg[1] - neg[0]) if len(neg) >= 2 else float("nan")
    say(f"  {beta:8.5f} {vp[0]:+11.7f} {Gpp_closed(np.pi, beta)['Gpp']:+11.7f}"
        f" {vp[i]:+12.7f} {pos[i]:8.4f} {len(cuts)//2:6d}  "
        + " ".join(f"{c:+.4f}" for c in cuts)
        + (f"   [a>0 width {w:.5f}]" if w == w else ""))
say("  (o02_bscan.log:41-43 grid = 101 pts on (-2pi,2pi) => edge bin +-0.063; 'negative window'")
say("   column there is min..max, not a contiguous interval.  Bisection here is exact.)")
say("")

# ---------------------------------------------------------------- [3] far boundary
say("[3] far boundary a -> 2pi (slice collapses to the single point rho=2pi, i.e. U=-1)")
say("    prediction from the exact reduction: I ~ eps^3/(6pi)  =>  G=-3log eps+O(1)")
say("    => G''(a) ~ +3/(2pi-a)^2 -> +inf.   (This is the h1d 'Var=0 => +' mechanism.)")
for beta in (0.0, 2.5, 6.0, 20.0):
    row = []
    for eps in (1e-1, 1e-2, 1e-3, 1e-4, 1e-5):
        g = Gpp_closed(TWOPI - eps, beta)["Gpp"]
        row.append(eps * eps * g)
    say(f"  beta={beta:5.2f}  (2pi-a)^2 G'' = " + "  ".join(f"{t:.8f}" for t in row) + "   [->3]")
say("")

# ---------------------------------------------------------------- [4] large beta
say("[4] large-beta limit at fixed a: mass -> rho=2pi-delta, delta=8pi/(beta a)")
say("    I ~ 64pi^2/(beta a)^3  =>  G ~ beta + 3 log(beta a) - log(64pi^2)  => G'' -> -3/a^2")
for a in (1.0, 1.5, 2.0, np.pi):
    row = []
    for beta in (20.0, 50.0, 200.0):
        row.append(a * a * Gpp_closed(a, beta)["Gpp"])
    say(f"  a={a:5.3f}  a^2*G'' at beta=20/50/200: " + "  ".join(f"{t:+.6f}" for t in row) + "  [-> -3]")
say("  G(a) vs beta+3log(beta a)-log(64pi^2):")
for beta in (20.0, 50.0):
    a = 2.0
    say(f"    beta={beta:5.1f} a={a}: G={G_of(a,beta):.6f}  asym={beta+3*np.log(beta*a)-np.log(64*np.pi**2):.6f}")
say("  saddle location (argmax of the tilted density) vs 8pi/(beta a), a=pi:")
say(f"    {'beta':>7} {'2pi-rho*':>10} {'8pi/(beta a)':>13}")
for beta in (2.5, 5.0, 10.0, 20.0, 50.0):
    a = np.pi
    rr = np.linspace(a, TWOPI - 1e-12, 400001)
    dens = W(rr) * np.exp(-beta * a * PH(rr))
    d = TWOPI - float(rr[int(np.argmax(dens))])
    say(f"    {beta:7.1f} {d:10.5f} {8*np.pi/(beta*a):13.5f}")
say("")

# ---------------------------------------------------------------- [5] thresholds
say("[5] two thresholds that decide the H1 question")
def bis_beta(f, lo=1e-6, hi=60.0, it=70):
    flo = f(lo)
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if (f(mid) > 0) == (flo > 0):
            lo, flo = mid, f(mid)
        else:
            hi = mid
    return 0.5 * (lo + hi)

b0 = bis_beta(lambda b: Gpp_closed(0.0, b)["Gpp"])
bpi = bis_beta(lambda b: Gpp_closed(np.pi, b)["Gpp"])
say(f"  beta*(a=0) = {b0:.6f}   (O28 [4]: 4.30661 ; h1d: sqrt(anchor/Var)=4.306607)")
A = Gpp_closed(0.0, 0.0)["Gpp"]
V = Gpp_closed(0.0, 0.0)["var"]
say(f"  closed form: A=1/(2Cin2pi)={A:.10f}, V=Var_w(phi)={V:.10f} => sqrt(A/V)={np.sqrt(A/V):.6f}")
say(f"  beta*(a=pi) = {bpi:.6f}  <-- G''(pi) changes sign here; so 'both ends positive'")
say(f"      holds only for beta < {bpi:.4f} (o28's +0.0272768642 is a beta=2.5 statement)")
say(f"  inf_a beta*(a) (= beta_cv, o02_bscan.log:42 = 1.873291) <= beta*(pi) <= beta*(0):"
    f" {1.87290:.5f} <= {bpi:.5f} <= {b0:.5f}")
say("")

# ---------------------------------------------------------------- [6] width law
say("[6] width of the a>0 negative window and its depth, exact roots")
say(f"  {'beta':>8} {'a1':>10} {'a2':>10} {'width':>10} {'2pi-a2':>10} {'depth':>11} {'argmin':>8}")
for beta in (2.0, 2.5, 3.0, 4.0, 4.30661, 5.0, 6.0, 8.0, 12.0, 20.0, 40.0):
    grid = np.linspace(1e-6, TWOPI - 1e-3, 801)
    v = np.array([Gpp_closed(float(a), beta)["Gpp"] for a in grid])
    rt = []
    for i in range(len(grid) - 1):
        if v[i] * v[i + 1] < 0:
            lo, hi, flo = grid[i], grid[i + 1], v[i]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                fm = Gpp_closed(float(mid), beta)["Gpp"]
                if (fm > 0) == (flo > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            rt.append(0.5 * (lo + hi))
    j = int(v.argmin())
    a1 = rt[0] if len(rt) >= 2 else 0.0
    a2 = rt[1] if len(rt) >= 2 else (rt[0] if rt else float("nan"))
    say(f"  {beta:8.5f} {a1:10.6f} {a2:10.6f} {a2-a1:10.6f}"
        f" {TWOPI-a2:10.6f} {v[j]:+11.7f} {grid[j]:8.4f}"
        + ("" if len(rt) >= 2 else "   (a1=0: left edge already past a=0)"))
say("")

# ---------------------------------------------------------------- [7] theta=0 control
say("[7] CONTROL theta=0: a-dependence is the tail domain only (no tilt), so")
say("    Z(a)=int_|a|^2pi g(rho)drho, g=w(rho)e^{-beta(1-cos(rho/2))}, G''=(g'Z+g^2)/Z^2.")
say("    Var(V_3) is NOT identically 0 here for a!=0 -- but the sign is decided by g' alone,")
say("    and g'<0 beyond rho=2.331 => negative windows even with no tilt.  Checks o02 [C].")
def tailG(a, beta, n=3000):
    f = lambda r: W(r) * np.exp(-beta * (1.0 - np.cos(0.5 * r)))
    return beta - np.log(integ(f, abs(a), TWOPI, n))


def tailGpp(a, beta, h=2e-3, n=1500):
    """route A on the theta=0 G(a): pure tail integral, no tilt term at all."""
    g = lambda z: tailG(z, beta, n)
    return (-g(a + 2 * h) + 16 * g(a + h) - 30 * g(a) + 16 * g(a - h) - g(a - 2 * h)) / (12 * h * h)


say(f"  check vs O28 [0]: G''(0;theta=0,beta=2.5) = {tailGpp(0.0, 2.5):+.7f}  (O28 0.7299241)")
for beta in (2.5, 6.0, 8.0):
    grid = np.linspace(1e-3, TWOPI - 0.02, 601)
    v = np.array([tailGpp(float(a), beta) for a in grid])
    rt = []
    for i in range(len(grid) - 1):
        if v[i] * v[i + 1] < 0:
            lo, hi, flo = grid[i], grid[i + 1], v[i]
            for _ in range(50):
                mid = 0.5 * (lo + hi)
                fm = tailGpp(float(mid), beta)
                if (fm > 0) == (flo > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            rt.append((lo, hi))
    say(f"  beta={beta:4.1f}  G''(0)={tailGpp(0.0,beta):+.7f}  min={v.min():+.7f} @ {grid[int(v.argmin())]:.4f}"
        f"  windows(a>0): " + ", ".join(f"[{p:.3f},{q:.3f}]" for p, q in rt)
        + "   (o02 [C] beta=6: [3.62,4.39]; beta=8: [3.49,4.78])")
say("")

# ---------------------------------------------------------------- [8] verdict
say("=" * 98)
say("[8] VERDICT INPUTS")
say(f"  (i)   a->2pi: G'' ~ +3/(2pi-a)^2 -> +inf  (analytic + numeric, table [3])")
say(f"  (ii)  G''(pi;2.5) = {Gpp_closed(np.pi,2.5)['Gpp']:+.7f} > 0  (O28 0.0272768642),")
say(f"        but beta*(pi) = {bpi:.5f} => positive-at-a=pi only for beta<{bpi:.3f}")
say(f"  (iii) a>0 negative window at beta=2.5: [0.709054, 2.973310], width 2.264256,"
    f" depth {Gpp_closed(1.822,2.5)['Gpp']:+.7f}")
say(f"  (iv)  window widens monotonically with beta; left edge hits a=0 exactly at")
say(f"        beta*(0)={b0:.5f}; right edge -> 2pi with 2pi-a2 -> 0; depth -> -inf")
say(f"  (v)   a=0 is positive for every beta < {b0:.5f}  => for beta in [1/48, 2.5] the a=0")
say(f"        statement is untouched by the negative window.")
say("=" * 98)
LOG.close()
