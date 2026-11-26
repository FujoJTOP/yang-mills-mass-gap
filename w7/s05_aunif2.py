"""S05-AUNIF part 2: refine a* (H1 threshold), a_rig (rigorous-bound threshold), c2, GL-order stability."""
import numpy as np
from scipy.optimize import brentq

LOG = open(r"D:\Dev\math\w7\s05_aunif2.log", "w", encoding="utf-8")


def say(*xs):
    line = " ".join(str(z) for z in xs)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


TWOPI = 2.0 * np.pi
BETA, TH = 2.5, np.pi / 2.0
CC, TT = np.cos(0.5 * TH), np.sin(0.5 * TH)
K0 = 0.5 * CC / TT
_G = {}
D1C = -TT / 8.0 - CC ** 2 / (8.0 * TT)
D3C = TT / 480.0
D1Q = TT / 320.0 + CC ** 2 / (192.0 * TT)
D3Q = -11.0 * TT / 107520.0


def gl(n):
    if n not in _G:
        xx, wg = np.polynomial.legendre.leggauss(n)
        _G[n] = (0.5 * (xx + 1.0), 0.5 * wg)
    return _G[n]


def dcoef(r):
    r = np.asarray(r, float)
    sm = r < 0.02
    rr = np.where(sm, 1.0, r)
    s, co = np.sin(0.5 * rr), np.cos(0.5 * rr)
    c1v = 3.0 * TT * (co / (2.0 * rr ** 2) - s / rr ** 3) - (CC ** 2 / (4.0 * TT)) * s / rr
    c3v = TT * (3.0 * s / rr ** 5 - 1.5 * co / rr ** 4 - s / (4.0 * rr ** 3))
    return np.where(sm, D1C + D1Q * r ** 2, c1v), np.where(sm, D3C + D3Q * r ** 2, c3v)


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


def stats(av, N):
    lo = abs(av)
    tt, wt = gl(N)
    rho = lo + (TWOPI - lo) * tt
    wt = (TWOPI - lo) * wt
    dens = wfun(rho) * np.exp(-BETA * (1.0 - CC * np.cos(0.5 * rho) + TT * av * wfun(rho)))
    p = wt * dens
    p = p / p.sum()
    Gv = BETA * ((CC / 2.0) * (av / rho) * np.sin(0.5 * rho)
                 + TT * (wfun(rho) + (av * av / rho) * wpf(rho)))
    d1v, d3v = dcoef(rho)
    Delt = av * d1v + av ** 3 * d3v
    EG = np.sum(p * Gv)
    return EG, np.sum(p * (Gv - EG) ** 2), np.sum(p * Delt)


say("=" * 100)
say("S05-AUNIF part2:  thresholds, c2, GL-order stability   (beta=2.5, theta=pi/2)")
say("=" * 100)

M1, M3 = 1.7677669530e-01, 1.4731391275e-03   # sup norms (part 1, 4e5-point scan)

say("")
say("-- C1  GL-order stability of a* (root of  gap(a) = [K0 E[G]-Var(G)] - beta|E[Delta]| )")
for N in (1500, 3000, 6000, 9000):
    def gap(av, N=N):
        EG, VG, ED = stats(av, N)
        return (K0 * EG - VG) - BETA * abs(ED)
    a = np.linspace(1.0, 1.6, 25)
    v = [gap(x) for x in a]
    st = None
    for i in range(len(a) - 1):
        if v[i] > 0 >= v[i + 1]:
            st = brentq(gap, a[i], a[i + 1], xtol=1e-13)
            break
    say(f"   n={N:5d}   a* = {st:.10f}   gap(1.2)={gap(1.2):+.6e}  gap(1.3)={gap(1.3):+.6e}")

say("")
say("-- C2  rigorous-bound threshold a_rig (root of  beta(aM1+a^3M3) = M(a) )")
for N in (1500, 6000):
    def gapr(av, N=N):
        EG, VG, ED = stats(av, N)
        return (K0 * EG - VG) - BETA * (av * M1 + av ** 3 * M3)
    a = np.linspace(0.8, 1.4, 25)
    v = [gapr(x) for x in a]
    st = None
    for i in range(len(a) - 1):
        if v[i] > 0 >= v[i + 1]:
            st = brentq(gapr, a[i], a[i + 1], xtol=1e-13)
            break
    say(f"   n={N:5d}   a_rig = {st:.10f}")

say("")
say("-- C3  c2 from  E[Delta]/a = c1 + c2 a + O(a^2)")
c1 = None
vals = []
for av in (1e-2, 5e-3, 2.5e-3, 1.25e-3, 6.25e-4):
    _, _, ED = stats(av, 4000)
    vals.append((av, ED / av))
    say(f"   a={av:.3e}   E[Delta]/a = {ED/av:+.12e}")
c1 = vals[-1][1]
say(f"   Richardson (a->0)  c1 ~ {c1:+.10e}")
c2s = [(vals[i][1] - vals[i + 1][1]) / (vals[i][0] - vals[i + 1][0]) for i in range(len(vals) - 1)]
say(f"   c2 estimates (finite differences) : " + "  ".join(f"{z:+.4e}" for z in c2s))
say(f"   => E[Delta] = {c1:+.6e}*a {c2s[-1]:+.4e}*a^2 + O(a^3)")

say("")
say("-- C4  table for the report, n=6000 : E[Delta], beta|E[Delta]|, M(a), gap, ratio")
say(f"  {'a':>5} {'E[Delta]':>13} {'beta|ED|':>12} {'M(a)':>11} {'gap':>11} {'ratio':>7}  H1")
for av in [0.0, 0.1, 0.25, 0.5, 0.75, 1.0, 1.2, 1.25, 1.28868, 1.3, 1.5]:
    EG, VG, ED = stats(av, 6000)
    M = K0 * EG - VG
    g = M - BETA * abs(ED)
    say(f"  {av:5.3f} {ED:+.6e} {BETA*abs(ED):.6e} {M:+.6e} {g:+.6e} "
        f"{BETA*abs(ED)/M:7.4f}  {'OK' if g > 0 else 'FAIL'}")

say("")
say("-- C5  is the whole H1 defect really the Delta term?  Var(G) <= K0 E[G] + beta E[Delta]")
say("     (E[H] = K0 E[G] + beta E[Delta] exact by r07 (4.1a); identical to B6)")
say("")
say("DONE")
LOG.close()
