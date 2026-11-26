"""S05-AUNIF:  E[Delta] bound (with sign) for  H = K0*G + beta*Delta ,  K0 = cot(theta/2)/2.

Delta = [ t*d_aa - c^2/(4t) ] (a*w) ,  w(rho) = sin(rho/2)/rho ,  rho in [|a|,2pi]
measure p ∝ w exp(-beta*S) ,  S = 1 - c cos(rho/2) + t a w      (r07 convention)

Part A: sympy exact -- Delta = a*delta1 + a^3*delta3 (EXACT, finite), delta1/delta3 closed
        forms, limits at rho->0 (no 1/rho^3 pole: cancellation), c1 = <delta1>_0 structure.
Part B: numerics -- E[Delta], beta*E[Delta], margin K0*E[G]-Var(G), crossing a*,
        sup-norm rigorous bound, sign scan.
"""
import numpy as np
import sympy as sp

LOG = open(r"D:\Dev\math\w7\s05_aunif.log", "w", encoding="utf-8")


def say(*xs):
    line = " ".join(str(z) for z in xs)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


say("=" * 104)
say("S05-AUNIF  H = K0 G + beta Delta ;  E[Delta] with sign ;  beta=2.5 theta=pi/2 a in [0,1.5]")
say("=" * 104)

# ============================ PART A : sympy ============================
rho, a = sp.symbols("rho a", positive=True)
c, t = sp.symbols("c t", positive=True)
w = sp.sin(rho / 2) / rho
wp = sp.diff(w, rho)
wpp = sp.diff(w, rho, 2)

Delta = sp.simplify(t * (3 * a * wp / rho - a ** 3 * wp / rho ** 3 + a ** 3 * wpp / rho ** 2)
                    - c ** 2 / (4 * t) * (a * w))
d1 = sp.simplify(3 * t * wp / rho - c ** 2 / (4 * t) * w)
d3 = sp.simplify(t * (wpp / rho ** 2 - wp / rho ** 3))

say("")
say("-- A1  EXACT expansion  Delta = a*delta1(rho) + a^3*delta3(rho)   (no a^2, no a^4, no remainder)")
say("   Delta - (a*d1 + a^3*d3)  simplified :", sp.simplify(Delta - (a * d1 + a ** 3 * d3)))
say("   delta1 =", sp.simplify(d1))
say("   delta3 =", sp.simplify(d3))
say("   delta1 over common denom: (-c^2*rho^2*sin(rho/2) + 6t^2[rho cos(rho/2) - 2 sin(rho/2)])/(4 t rho^3)")
say("   delta3 over common denom: t[-rho^2 sin(rho/2) - 6 rho cos(rho/2) + 12 sin(rho/2)]/(4 rho^5)")
say("   [series rho->0] delta1 =", sp.series(d1, rho, 0, 4))
say("   [series rho->0] delta3 =", sp.series(d3, rho, 0, 6))
say("   limit rho->0 delta1 =", sp.limit(d1, rho, 0), "   ( = -t/8 - c^2/(8t) )")
say("   limit rho->0 delta3 =", sp.limit(d3, rho, 0), "   ( = t/480 )")
say("   => w-singularity at u=0 does NOT produce a pole: 1/rho^3 and 1/rho^4 terms cancel")
say("   [parity] d1(-rho)-d1(rho) =", sp.simplify(d1.subs(rho, -rho) - d1),
    " ; d3(-rho)-d3(rho) =", sp.simplify(d3.subs(rho, -rho) - d3))
say("   [d/dt-independence] Delta/Delta at a->0 slope: delta1 = d(Delta)/da |_a=0")
say("   check: d(Delta)/da at a=0 - d1 =", sp.simplify(sp.diff(Delta, a).subs(a, 0) - d1))

# u-form
U = sp.symbols("u", positive=True)
gfun = sp.sin(U) / U
say("")
say("-- A2  u-form (u=rho/2, w=g/2, w'=g'/4, w''=g''/8)")
say("   delta1 = (3t/8)*g'(u)/u - (c^2/(8t))*g(u)   check :",
    sp.simplify(d1.subs(rho, 2 * U) - (sp.Rational(3, 8) * t * sp.diff(gfun, U) / U
                                        - c ** 2 / (8 * t) * gfun)))
say("   delta3 = t*[ g''(u)/(32u^2) - g'(u)/(32u^3) ]  check :",
    sp.simplify(d3.subs(rho, 2 * U) - t * (sp.diff(gfun, U, 2) / (32 * U ** 2)
                                          - sp.diff(gfun, U) / (32 * U ** 3))))

# small-a limits of the measure pieces (for c1 / c2 bookkeeping)
say("")
say("-- A3  why c1 = <delta1>_0  (a->0 bookkeeping)")
say("   p = w e^{-beta S}/Z ,  S = S0 + t a w ,  S0 = 1 - c cos(rho/2) ,  rho in [a,2pi] (a>=0)")
say("   Z(a)  = Z0 - a*F0 - beta t a Z0 <w>_0 + O(a^2) ,  F0 = w(0)e^{-beta S0(0)} = (1/2)e^{-beta(1-c)}")
say("   N(a)  = int_a^2pi (a d1 + a^3 d3) F e^{-beta t a w} d rho")
say("         = a D1 - a^2[ beta t int d1 F w + (d1 F)(0) ] + a^3 D3 + O(a^4)")
say("   => E[Delta] = c1 a + c2 a^2 + ... ,  c1 = D1/Z0 = <delta1>_0  (EXACT leading coefficient)")
say("      the -a^2 (d1F)(0) piece is the moving cutoff rho>=|a| ; it is O(a^2), NOT O(a) : no c1 contamination")

# ============================ PART B : numerics ============================
TWOPI = 2.0 * np.pi
BETA, TH = 2.5, np.pi / 2.0
CC, TT = np.cos(0.5 * TH), np.sin(0.5 * TH)
K0 = 0.5 * np.cos(0.5 * TH) / np.sin(0.5 * TH)
_G = {}


def gl(n):
    if n not in _G:
        xx, wg = np.polynomial.legendre.leggauss(n)
        _G[n] = (0.5 * (xx + 1.0), 0.5 * wg)
    return _G[n]


D1C = -TT / 8.0 - CC ** 2 / (8.0 * TT)   # delta1(0)
D3C = TT / 480.0                         # delta3(0)
D1Q = TT / 320.0 + CC ** 2 / (192.0 * TT)
D3Q = -11.0 * TT / 107520.0


def dcoef(r):
    """delta1, delta3 with a Taylor branch near rho=0 (both are analytic EVEN functions)."""
    r = np.asarray(r, float)
    sm = r < 0.02
    rr = np.where(sm, 1.0, r)
    s, co = np.sin(0.5 * rr), np.cos(0.5 * rr)
    c1v = 3.0 * TT * (co / (2.0 * rr ** 2) - s / rr ** 3) - (CC ** 2 / (4.0 * TT)) * s / rr
    c3v = TT * (3.0 * s / rr ** 5 - 1.5 * co / rr ** 4 - s / (4.0 * rr ** 3))
    out1 = np.where(sm, D1C + D1Q * r ** 2, c1v)
    out3 = np.where(sm, D3C + D3Q * r ** 2, c3v)
    return out1, out3


def wfun(r):
    r = np.asarray(r, float)
    out = np.empty_like(r)
    m = r != 0.0
    out[m] = np.sin(0.5 * r[m]) / r[m]
    out[~m] = 0.5
    return out


def wpf(r):
    rr = np.where(r == 0.0, 1.0, r)
    return (0.5 * rr * np.cos(0.5 * rr) - np.sin(0.5 * rr)) / rr ** 2


def grid(av, N):
    lo = abs(av)
    tt, wt = gl(N)
    rho = lo + (TWOPI - lo) * tt
    wt = (TWOPI - lo) * wt
    dens = wfun(rho) * np.exp(-BETA * (1.0 - CC * np.cos(0.5 * rho) + TT * av * wfun(rho)))
    p = wt * dens
    return rho, p / p.sum()


def fields(av, rho):
    w_ = wfun(rho)
    wp_ = wpf(rho)
    Gv = BETA * ((CC / 2.0) * (av / rho) * np.sin(0.5 * rho) + TT * (w_ + (av * av / rho) * wp_))
    d1v, d3v = dcoef(rho)
    Delt = av * d1v + av ** 3 * d3v
    return Gv, Delt


def series_c(a_vals, N=6000):
    out = []
    for av in a_vals:
        rho, p = grid(av, N)
        Gv, Delt = fields(av, rho)
        EG, VG, ED = np.sum(p * Gv), np.sum(p * (Gv - np.sum(p * Gv)) ** 2), np.sum(p * Delt)
        out.append((av, EG, VG, ED, ED / av if av else 0.0, K0 * EG - VG))
    return out


say("")
say("=" * 104)
say("PART B  measure p ∝ w exp(-beta S), GL n=6000  (r07 convention, rho in [|a|,2pi])")
say("=" * 104)

# ---- B1: reproducibility of r07's E[G], Var(G) + E[Delta] and the margin ----
say("")
say("-- B1  E[G], Var(G), E[Delta], beta*E[Delta], margin M(a)=K0*E[G]-Var(G), ratio")
say("    (r07 ref: E[G]=0.804775/0.897072/0.988335/1.018545 ; Var(G)=1.2215e-2/1.8744e-2/"
    "3.1811e-2/4.8047e-2  @ a=0/0.3/0.75/1.2)")
say(f"  {'a':>6} {'E[G]':>11} {'Var(G)':>11} {'E[Delta]':>12} {'beta*E[Delta]':>13} "
    f"{'beta|E[Delta]|':>14} {'M(a)':>10} {'ratio':>7}")
rows = []
for av in [0.0, 0.3, 0.75, 1.2]:
    rho, p = grid(av, 6000)
    Gv, Delt = fields(av, rho)
    EG = np.sum(p * Gv)
    VG = np.sum(p * (Gv - EG) ** 2)
    ED = np.sum(p * Delt)
    M = K0 * EG - VG
    rows.append((av, EG, VG, ED, M))
    say(f"  {av:6.3f} {EG:+.6e} {VG:.6e} {ED:+.6e} {BETA*ED:+.6e} {BETA*abs(ED):.6e} {M:+.6e} "
        f"{BETA*abs(ED)/M:7.3f}")

# ---- B2: c1 -- direct limit, and the a=0 average <delta1>_0 ----
say("")
say("-- B2  c1 = lim_{a->0} E[Delta]/a   vs   <delta1>_0  (a=0 measure, direct quadrature)")
rho0, p0 = grid(0.0, 6000)
d10, d30 = dcoef(rho0)
c1_direct = np.sum(p0 * d10)
say(f"   <delta1>_0                      = {c1_direct:+.10e}")
say(f"   <delta3>_0                      = {np.sum(p0*d30):+.10e}")
say(f"   delta1(0) = -t/8 - c^2/(8t)     = {D1C:+.10e}")
say(f"   delta3(0) =  t/480              = {D3C:+.10e}")
say(f"   min delta1 on (0,2pi]           = {np.min(d10):+.6e}   max = {np.max(d10):+.6e}  "
    f"(sign-definite? {'YES (<=0)' if np.max(d10) <= 1e-14 else 'NO'})")
say(f"   min delta3 on (0,2pi]           = {np.min(d30):+.6e}   max = {np.max(d30):+.6e}")
for av in [1e-1, 1e-2, 1e-3, 1e-4]:
    _, _, _, ED, rat, _ = series_c([av])[0]
    say(f"   a={av:<7.0e}  E[Delta]={ED:+.8e}   E[Delta]/a={rat:+.10e}   "
        f"ratio to c1 = {rat/c1_direct:.8f}   (c2*a correction ~ a)")

# ---- B3: dense a-grid: margin vs beta|E[Delta]| ; crossing a* ----
say("")
say("-- B3  dense grid a in [0,1.5] : margin vs beta|E[Delta]|")
say(f"  {'a':>6} {'E[Delta]':>13} {'beta|E[Delta]|':>15} {'M(a)':>11} {'M-beta|ED|':>12} {'ratio':>8}")
grid_a = np.linspace(0.0, 1.5, 31)
tab = []
for av in grid_a:
    rho, p = grid(av, 4000)
    Gv, Delt = fields(av, rho)
    EG = np.sum(p * Gv)
    VG = np.sum(p * (Gv - EG) ** 2)
    ED = np.sum(p * Delt)
    M = K0 * EG - VG
    tab.append((av, ED, M))
    if av <= 1.5:
        say(f"  {av:6.3f} {ED:+.6e} {BETA*abs(ED):.6e} {M:+.6e} {M-BETA*abs(ED):+.6e} "
            f"{BETA*abs(ED)/M:8.4f}")

from scipy.optimize import brentq  # noqa: E402


def gap(av):
    rho, p = grid(av, 4000)
    Gv, Delt = fields(av, rho)
    EG = np.sum(p * Gv)
    VG = np.sum(p * (Gv - EG) ** 2)
    ED = np.sum(p * Delt)
    return (K0 * EG - VG) - BETA * abs(ED)


ascan = np.linspace(1e-6, 3.0, 61)
vals = [gap(x) for x in ascan]
say("")
say(f"   gap(a) at a=1.5 : {gap(1.5):+.6e}    (positive => H1 holds at a=1.5)")
star = None
for i in range(len(ascan) - 1):
    if vals[i] > 0 and vals[i + 1] <= 0:
        star = brentq(gap, ascan[i], ascan[i + 1], xtol=1e-12)
        say(f"   *** crossing (margin = beta|E[Delta]|) at a* = {star:.8f} ***")
        break
if star is None:
    say(f"   no crossing on [1e-6,3.0] : min gap = {min(vals):+.6e} at a={ascan[int(np.argmin(vals))]:.4f}")

# ---- B4: rigorous sup-norm bound  |E[Delta]| <= a*||d1||_inf + a^3*||d3||_inf ----
say("")
say("-- B4  RIGOROUS a-uniform bound (no asymptotics):  |E[Delta]| <= a*M1 + a^3*M3")
rr = np.concatenate([np.linspace(1e-12, 1e-2, 2000), np.linspace(1e-2, TWOPI, 400000)])
b1, b3 = dcoef(rr)
M1, M3 = np.max(np.abs(b1)), np.max(np.abs(b3))
say(f"   ||delta1||_inf = {M1:.10e} at rho = {rr[int(np.argmax(np.abs(b1)))]:.6f}")
say(f"   ||delta3||_inf = {M3:.10e} at rho = {rr[int(np.argmax(np.abs(b3)))]:.6f}")
say(f"   {'a':>6} {'beta(aM1+a^3M3)':>18} {'M(a)':>11} {'slack':>12}  holds?")
for av in np.linspace(0.05, 1.5, 30):
    rho, p = grid(av, 4000)
    Gv, _ = fields(av, rho)
    EG = np.sum(p * Gv)
    VG = np.sum(p * (Gv - EG) ** 2)
    M = K0 * EG - VG
    bd = BETA * (av * M1 + av ** 3 * M3)
    say(f"   {av:6.3f} {bd:.6e} {M:+.6e} {M-bd:+.6e}  {'YES' if bd < M else 'NO'}")

# ---- B5: sign scan of E[Delta] (does it ever flip?) ----
say("")
say("-- B5  sign scan of E[Delta] on a in (0,3]  (does E[Delta] change sign?)")
neg = 0
firstpos = None
for av in np.linspace(0.01, 3.0, 60):
    rho, p = grid(av, 3000)
    _, Delt = fields(av, rho)
    ED = np.sum(p * Delt)
    if ED < 0:
        neg += 1
    elif firstpos is None:
        firstpos = av
say(f"   #negative = {neg} / 60 ;  first positive a = {firstpos}")
say(f"   => E[Delta] is {'NEGATIVE on (0,3]' if firstpos is None else 'sign-changing'}")

# ---- B6: the criterion itself (H1): Var(G) <= K0 E[G] + beta E[Delta] ----
say("")
say("-- B6  H1 itself:  Var(G) <= K0 E[G] + beta E[Delta]  <=>  beta|E[Delta]| <= M(a)")
say(f"  {'a':>6} {'LHS Var(G)':>13} {'RHS K0EG+bED':>14} {'RHS-margin form':>16} {'H1':>5}")
for av in np.linspace(0.0, 1.5, 16):
    rho, p = grid(av, 4000)
    Gv, Delt = fields(av, rho)
    EG = np.sum(p * Gv)
    VG = np.sum(p * (Gv - EG) ** 2)
    ED = np.sum(p * Delt)
    rhs = K0 * EG + BETA * ED
    say(f"  {av:6.3f} {VG:.6e} {rhs:+.6e} {rhs-VG:+.6e} {'OK' if VG <= rhs else 'FAIL':>5}")

say("")
say("DONE")
LOG.close()
