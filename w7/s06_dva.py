"""S06-DA -- a-dependence of S(a)=E[H+Phi''|a], D(a)=Var(G+Phi'|a), gap(a)=S-D=G''(a).

Model (verbatim H1A-EXACT 2.1 = O28/O09/q11), single plaquette SU(2), beta=2.5:
  S_W = beta(1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) xi3/rho), J = sin^2(rho/2)/rho^2
  slice xi3 = a, |xi| < 2pi => rho in [|a|, 2pi]  (axisymmetry: the slice integral is 1D)
  measure  W(rho) drho = 4 sin^2(rho/2)/rho drho   (= 4 * rho * J * drho ; the 4 is q08's Phi convention)
  S = E[H + Phi''|a],  D = Var(G + Phi'|a),  gap = S - D = G''(a) = E[V33] - Var(V3)

Routes:  B = analytic fixed-s coordinate derivatives (O28 formulas, cancellation-free Phi'')
         A = 5-point stencil in a on G(a) = -log Z(a)  (touches neither H nor Phi -> checks gap=S-D)

Outputs: [0] validation vs q11/r03/o28/q08, [1] the a-scan for th=pi and pi/2,
         [2] monotonicity, [3] interior extrema (golden section), [4] a=0 local-max test
         (including the a<0 side), [5] convergence.
"""
import sys
import numpy as np

LOG = open(r"D:\Dev\math\w7\s06_dva.log", "w", encoding="utf-8")


def say(*x):
    line = " ".join(str(t) for t in x)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


TWOPI = 2.0 * np.pi
G = {}


def gl(n):
    if n not in G:
        x, w = np.polynomial.legendre.leggauss(n)
        G[n] = (0.5 * (x + 1.0), 0.5 * w)
    return G[n]


def dens(a, beta, th, n=3000):
    """nodes + quadrature weights of the slice measure, plus every a-derivative pointwise."""
    c, s = np.cos(th / 2.0), np.sin(th / 2.0)
    lo, hi = abs(a), TWOPI
    t, wt = gl(n)
    rho = lo + (hi - lo) * t
    wt = (hi - lo) * wt
    h = 0.5 * rho
    sh, ch = np.sin(h), np.cos(h)
    u = a / rho
    s2 = rho * rho - a * a                 # transverse radius^2
    r3 = s2 / rho ** 3                     # (1-u^2)/rho
    w = 4.0 * sh * sh / rho * np.exp(-beta * (1.0 - c * ch + s * sh * u)) * wt
    # action: H = d^2_{xi3} S_W at fixed transverse radius
    S3 = 0.5 * beta * c * sh * u + 0.5 * beta * s * ch * u * u + beta * s * sh * r3
    S33 = ((0.25 * beta * c * ch - 0.25 * beta * s * sh * u) * u * u
           + beta * s * ch * u * r3
           + (0.5 * beta * c * sh + 0.5 * beta * s * ch * u) * r3
           - 3.0 * beta * s * sh * u * s2 / rho ** 4)
    # jacobian: Phi = -log J ; Php = 1/h - cot h ; Phpp cancellation-free
    Php = 1.0 / h - ch / sh
    Phpp = (h * h / (sh * sh) - 1.0) / (2.0 * h * h)
    P3 = Php * u
    P33 = Phpp * u * u + Php * r3
    return w, S3, S33, P3, P33


def stats(a, beta, th, n=3000):
    w, S3, S33, P3, P33 = dens(a, beta, th, n)
    Z = w.sum()
    EH = (w * S33).sum() / Z
    EP2 = (w * P33).sum() / Z
    V3 = S3 + P3
    m = (w * V3).sum() / Z
    D = (w * (V3 - m) ** 2).sum() / Z          # two-pass: no variance cancellation
    varG = (w * (S3 - (w * S3).sum() / Z) ** 2).sum() / Z
    varP = (w * (P3 - (w * P3).sum() / Z) ** 2).sum() / Z
    cov = (w * (S3 - (w * S3).sum() / Z) * (P3 - (w * P3).sum() / Z)).sum() / Z
    S = EH + EP2
    return dict(Z=Z, EH=EH, EP2=EP2, S=S, D=D, gap=S - D, varG=varG, varP=varP, cov=cov,
                gapH=EH - varG - cov)


def G_of(a, beta, th, n=3000):
    return beta - np.log(stats(a, beta, th, n)["Z"])


def Gpp_st(a, beta, th, hstep=1e-3, n=1500):
    g = lambda z: G_of(z, beta, th, n)
    return (-g(a + 2 * hstep) + 16 * g(a + hstep) - 30 * g(a) + 16 * g(a - hstep) - g(a - 2 * hstep)) / (12 * hstep ** 2)


def golden(f, lo, hi, it=200):
    p = (np.sqrt(5.0) - 1.0) / 2.0
    a1, a2 = lo + (1 - p) * (hi - lo), lo + p * (hi - lo)
    f1, f2 = f(a1), f(a2)
    for _ in range(it):
        if f1 < f2:
            hi, a2, f2 = a2, a1, f1
            a1 = lo + (1 - p) * (hi - lo)
            f1 = f(a1)
        else:
            lo, a1, f1 = a1, a2, f2
            a2 = lo + p * (hi - lo)
            f2 = f(a2)
    x = 0.5 * (lo + hi)
    return x, f(x)


BETA = 2.5
say("=" * 104)
say("S06-DA  S(a)=E[H+Phi''|a], D(a)=Var(G+Phi'|a), gap=S-D=G''(a)   beta=2.5, single plaquette SU(2)")
say(f"numpy {np.__version__}  python {sys.version.split()[0]}  GL n=3000 (convergence in [5])")
say("=" * 104)
say("")

# ---------------------------------------------------------------- [0] validation
say("[0] VALIDATION (a=0: q11/r03 exact intervals; a=pi, 1.9: q08=o28 published values)")
for th, want in ((np.pi, (0.20511529712307314798, 0.069120586081084936946)),
                 (np.pi / 2, (0.55818364272458720550, 0.015239720370897441829))):
    r = stats(0.0, BETA, th)
    say(f"  th={th:.6f} a=0: S={r['S']:+.15f} (ref {want[0]:+.15f}, d={r['S']-want[0]:+.1e})"
        f"  D={r['D']:+.15f} (ref {want[1]:+.15f}, d={r['D']-want[1]:+.1e})"
        f"  gap={r['gap']:+.15f}")
ref = ((1.9, -0.0788201414), (np.pi, 0.0272768642), (0.709054, -0.0000008626), (2.973310, -0.0002036677))
say(f"  {'a':>9} {'gap(B)':>16} {'gap(A) stencil':>16} {'A-B':>10} {'o28/q08':>14} {'diff':>10}")
for a, want in ref:
    rb = stats(a, BETA, np.pi)
    ra = Gpp_st(a, BETA, np.pi, 1e-3, 2400)
    say(f"  {a:9.6f} {rb['gap']:+16.10f} {ra:+16.10f} {ra-rb['gap']:+10.1e} {want:+14.10f} {rb['gap']-want:+.1e}")
x, v = golden(lambda z: stats(z, BETA, np.pi)["gap"], 1.0, 2.5)
say(f"  min on [1,2.5]: gap={v:+.12f} at a={x:.8f}   (o28/o08 grid: -0.0791945935 @ 1.82212)")
say("")

# ---------------------------------------------------------------- [1] the scan
AG = np.linspace(0.0, 2.5, 251)
TAB = {}
for th, tag in ((np.pi, "pi"), (np.pi / 2.0, "pi/2")):
    rows = [stats(float(a), BETA, th) for a in AG]
    TAB[tag] = rows
    say(f"[1] SCAN  beta=2.5  theta={tag}   a in [0,2.5], {len(AG)} points")
    say(f"  {'a':>6} {'S=E[H+Phi2]':>13} {'D=Var(G+Phi1)':>14} {'gap=S-D':>12} {'E[H]':>10} {'E[Phi2]':>10}"
        f" {'VarG':>10} {'VarPhi1':>10} {'Cov':>10}")
    for i in range(0, len(AG), 5):
        r = rows[i]
        say(f"  {AG[i]:6.3f} {r['S']:13.8f} {r['D']:14.8f} {r['gap']:+12.8f} {r['EH']:+10.6f}"
            f" {r['EP2']:10.6f} {r['varG']:10.7f} {r['varP']:10.7f} {r['cov']:+10.7f}")
    say(f"  (full 251-point CSV at the end of this log)")

    # ------------------------------------------------------------ [2] monotonicity
    S = np.array([r["S"] for r in rows])
    D = np.array([r["D"] for r in rows])
    gap = np.array([r["gap"] for r in rows])
    say(f"[2] MONOTONICITY on [0,2.5], theta={tag}")
    for nm, v in (("S", S), ("D", D), ("gap", gap)):
        d = np.diff(v)
        neg, pos = int(np.sum(d < 1e-14)), int(np.sum(d > -1e-14))
        flips = int(np.sum(np.diff(np.sign(np.round(d, 13))) != 0))
        sign = "increase" if neg == 0 else ("decrease" if pos == 0 else "NOT monotone")
        say(f"  {nm:4s}: {sign:13s}  #d<0={neg:3d} #d>0={pos:3d}  sign flips in d={flips}"
            f"   min={v.min():+.8f}@{AG[v.argmin()]:.4f}  max={v.max():+.8f}@{AG[v.argmax()]:.4f}")

    # ------------------------------------------------------------ [3] interior extrema
    say(f"[3] INTERIOR EXTREMA of gap, theta={tag} (golden section on the analytic route)")
    crit = []
    ds = np.diff(gap) * np.diff(gap)
    for i in range(1, len(AG) - 1):
        if (gap[i] - gap[i - 1]) * (gap[i + 1] - gap[i]) < 0:
            crit.append((AG[i - 1], AG[i + 1]))
    for lo, hi in crit:
        xm, vm = golden(lambda z: stats(z, BETA, th)["gap"], lo, hi)
        xM, vM = golden(lambda z: -stats(z, BETA, th)["gap"], lo, hi)
        f = lambda z: stats(z, BETA, th)["gap"]
        which = "min" if vm < vM else "max"
        x, v = (xm, vm) if which == "min" else (xM, -vM)
        say(f"  extremum in [{lo:.4f},{hi:.4f}]: {which} gap={v:+.12f} at a={x:.8f}")
    g0 = stats(0.0, BETA, th)["gap"]
    say(f"  gap(0)={g0:+.12f}   gap(2.5)={gap[-1]:+.12f}")
    say("")

# ---------------------------------------------------------------- [4] is a=0 a local max?
say("[4] IS a=0 A LOCAL MAX?  (needs the a<0 side: the model is NOT even -- S depends on u=a/rho)")
for th, tag in ((np.pi, "pi"), (np.pi / 2.0, "pi/2")):
    AN = np.linspace(-0.5, 0.0, 26)
    gn = np.array([stats(float(a), BETA, th)["gap"] for a in AN])
    gp = np.array([stats(float(a), BETA, th)["gap"] for a in np.linspace(0.0, 0.5, 26)])
    say(f"  theta={tag}: gap(-0.5..0) = " + " ".join(f"{v:+.5f}" for v in gn[::5])
        + f"   |  gap(0..0.5) = " + " ".join(f"{v:+.5f}" for v in gp[::5]))
    dl = (gn[-1] - gn[-2]) / (AN[-1] - AN[-2])
    dr = (gp[1] - gp[0]) / 0.02
    say(f"    one-sided slopes at a=0:  gap'(0-)={dl:+.6f}   gap'(0+)={dr:+.6f}"
        f"   => {'local max' if dl > 0 > dr else ('left-endpoint max' if dr < 0 else 'not a max')}")
    say(f"    max on [-0.5,0.5] = {max(gn.max(), gp.max()):+.8f} at a="
        f"{np.concatenate([AN, gp*0 + AN])[np.argmax(np.concatenate([gn, gp]))]:.4f}")
    say("")

# ---------------------------------------------------------------- [4b] turning points, beyond-2.5
say("[4b] TURNING POINTS of S and D (theta=pi) and the rest of the circle for theta=pi/2")
for nm, key in (("S", "S"), ("D", "D")):
    xm, vm = golden(lambda z: stats(z, BETA, np.pi)[key], 1.0, 2.5)
    say(f"  theta=pi  {nm}(a): min={vm:+.12f} at a={xm:.8f}")
for th, tag in ((np.pi, "pi"), (np.pi / 2.0, "pi/2")):
    AX = np.linspace(2.5, TWOPI - 0.02, 120)
    vx = np.array([stats(float(a), BETA, th)["gap"] for a in AX])
    j = int(vx.argmin())
    say(f"  theta={tag:4s} gap on [2.5,2pi-0.02]: min={vx[j]:+.10f} at a={AX[j]:.5f}"
        f"  (a=2.5: {stats(2.5, BETA, th)['gap']:+.10f}, a=6.26: {vx[-1]:+.6f})")
    if vx[j] < stats(2.5, BETA, th)["gap"]:
        xm, vm = golden(lambda z: stats(z, BETA, th)["gap"], AX[max(j - 2, 0)], AX[min(j + 2, len(AX) - 1)])
        say(f"           -> interior min far side: gap={vm:+.10f} at a={xm:.6f}")
xp = np.linspace(1e-3, 2.0, 401)
say("  gap'(0+) by central difference on [-0.2,0.2] (smoothness through a=0, both theta):")
for th, tag in ((np.pi, "pi"), (np.pi / 2.0, "pi/2")):
    for h in (1e-2, 1e-3):
        d = (stats(h, BETA, th)["gap"] - stats(-h, BETA, th)["gap"]) / (2 * h)
        say(f"    theta={tag:4s} h={h:5.0e}: central gap'(0)={d:+.8f}")
say("")

# ---------------------------------------------------------------- [5] convergence + CSV
say("[5] CONVERGENCE (theta=pi): S/D/gap for n = 1500/3000/6000")
for a in (0.0, 0.709054, 1.82212, 2.5):
    rr = [stats(a, BETA, np.pi, n) for n in (1500, 3000, 6000)]
    say(f"  a={a:8.5f}  S=" + " ".join(f"{r['S']:+.12f}" for r in rr)
        + "  D=" + " ".join(f"{r['D']:+.12f}" for r in rr)
        + f"  gap drift={rr[2]['gap']-rr[1]['gap']:+.1e}")
say("")
say("CSV  theta,a,S,D,gap,E[H],E[Phi'']")
for th, tag in ((np.pi, "pi"), (np.pi / 2.0, "pi/2")):
    for a, r in zip(AG, TAB[tag]):
        say(f"CSV,{tag},{a:.6f},{r['S']:.12f},{r['D']:.12f},{r['gap']:.12f},{r['EH']:.12f},{r['EP2']:.12f}")
say("=" * 104)
LOG.close()
