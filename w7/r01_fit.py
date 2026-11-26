"""R01-HG part 2: what IS the a!=0 relation?  weighted (measure) least squares of
   H  against  G  and small deterministic bases, on the slice measure of Q04 (1.6).
Bases tested (all deterministic in rho):
   {} : H = K G                (pure proportionality -- q04's shape)
   {1}         : H = K G + C
   {w}         : H = K G + C w
   {1,w}       : H = K G + C0 + C1 w
Numbers reported: K, C's, weighted R^2 and weighted relative residual std.
"""
import numpy as np

LOG = open(r"D:\Dev\math\w7\r01_fit.log", "w", encoding="utf-8")


def say(*x):
    line = " ".join(str(t) for t in x)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


TWOPI = 2.0 * np.pi
_G = {}


def gl(n):
    if n not in _G:
        x, wg = np.polynomial.legendre.leggauss(n)
        _G[n] = (0.5 * (x + 1.0), 0.5 * wg)
    return _G[n]


def wfun(r):
    r = np.asarray(r, float)
    out = np.empty_like(r)
    m = r != 0.0
    out[m] = np.sin(0.5 * r[m]) / r[m]
    out[~m] = 0.5
    return out


def wpf(r):
    return (0.5 * r * np.cos(0.5 * r) - np.sin(0.5 * r)) / np.where(r == 0.0, 1.0, r) ** 2


def wppf(r):
    rr = np.where(r == 0.0, 1.0, r)
    return -np.sin(0.5 * r) / (4.0 * rr) - np.cos(0.5 * r) / rr ** 2 + 2.0 * np.sin(0.5 * r) / rr ** 3


def fields(a, rho, beta, th):
    c, t = np.cos(0.5 * th), np.sin(0.5 * th)
    u = a / rho
    w, wp, wpp = wfun(rho), wpf(rho), wppf(rho)
    G = beta * ((c / 2.0) * u * np.sin(0.5 * rho) + t * (w + (a * a / rho) * wp))
    H = beta * ((c / 2.0) * (np.sin(0.5 * rho) / rho + (a * a / (2.0 * rho ** 2)) * np.cos(0.5 * rho)
                             - (a * a / rho ** 3) * np.sin(0.5 * rho))
                + t * ((3.0 * a / rho) * wp - (a ** 3 / rho ** 3) * wp + (a ** 3 / rho ** 2) * wpp))
    return G, H


BETA, TH = 2.5, np.pi / 2.0
K0 = 0.5 / np.tan(TH / 2.0)
N = 6000
say("=" * 100)
say("R01-HG part 2 : weighted least squares, beta=2.5 theta=pi/2, measure w exp(-beta S) drho")
say("=" * 100)
for av in [0.0, 0.3, 0.75, 1.2]:
    lo = abs(av)
    tt, wt = gl(N)
    rho = lo + (TWOPI - lo) * tt
    wt = (TWOPI - lo) * wt
    G, H = fields(av, rho, BETA, TH)
    dens = wfun(rho) * np.exp(-BETA * (1.0 - np.cos(0.5 * TH) * np.cos(0.5 * rho)
                                       + np.sin(0.5 * TH) * av * wfun(rho)))
    p = wt * dens
    p = p / p.sum()
    w = wfun(rho)
    say("")
    say(f"a = {av}   rho in [{lo:.3f}, 2pi]   (weight mass on grid)")
    for name, basis in [("H = K G            ", []),
                        ("H = K G + C        ", [np.ones_like(rho)]),
                        ("H = K G + C w      ", [w]),
                        ("H = K G + C0 + C1 w", [np.ones_like(rho), w])]:
        cols = [G] + basis
        X = np.stack(cols, axis=1)
        Xw = X * np.sqrt(p)[:, None]
        yw = H * np.sqrt(p)
        coef, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        res = H - X @ coef
        r2 = 1.0 - float(np.sum(p * res ** 2) / np.sum(p * (H - np.sum(p * H)) ** 2))
        relstd = float(np.sqrt(np.sum(p * res ** 2) / np.sum(p * H ** 2)))
        txt = "  ".join(f"C{i}={c:+.8f}" for i, c in enumerate(coef[1:]))
        say(f"  {name}  K={coef[0]:+.8f}  {txt:36s}  R2={r2:.10f}  relres={relstd:.3e}")
    say(f"  [reference] K0=(1/2)cot(th/2)={K0:+.8f} ; E[H]/E[G]={np.sum(p*H)/np.sum(p*G):+.8f}"
        f" ; Cov/Var(G)={ (np.sum(p*G*H)-np.sum(p*G)*np.sum(p*H))/(np.sum(p*G*G)-np.sum(p*G)**2):+.8f}")
LOG.close()
