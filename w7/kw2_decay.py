"""KW2: which decay profiles survive  |W|<=1  +  CM(Bernstein)  ?

Checks (numerical; no proof is claimed from them):
 (1) CM scan of W(R)=exp(-c R^alpha) by Cauchy derivatives n<=8, alpha grid.
 (2) independent witness: alpha=1/2 exact Levy density, quadrature vs exp(-sqrt R).
 (3) engine of the verdict: x^alpha = (alpha/Gamma(1-alpha)) Int (1-e^-xt) t^-1-alpha dt.
 (4) alpha>1 failure: explicit sign of W'' and the threshold R*.
 (5) log-power, EXACT profile: realness defect (branch) + endpoint limits R->0.
 (6) log-power, ASYMPTOTIC profile: rho(t)=t^-1 exp(-c|ln t|^beta), Laplace W in
     log domain, check the ln-power law and its exponent.
 (7) R^-p (ln R)^q: monotonicity for q>0, CM scan for q<0.
 (8) allowed set is big: gapped and gapless members side by side.
 (9) gap criterion: inf supp rho = rate.
Only prints; the log is the shell redirect.
"""
import math
import numpy as np

TH = 2 * np.pi * np.arange(192) / 192


def cauchy(f, n, R, rel=0.3):
    """f^(n)(R) by Cauchy integral on |z-R| = rel*R (spectral accuracy)."""
    r = rel * R
    z = R + r * np.exp(1j * TH)
    return (math.factorial(n) / r ** n) * np.mean(np.asarray(f(z), dtype=complex) * np.exp(-1j * n * TH))


def cm_scan(f, Rg, nmax=8, rel=0.3):
    """signed CM test: (-1)^n Re f^(n) >= 0; also returns realness defect of f^(n)."""
    fails, realdef = [], 0.0
    for n in range(1, nmax + 1):
        for R in Rg:
            v = cauchy(f, n, R, rel)
            scale = math.factorial(n) / (rel * R) ** n * abs(complex(f(complex(R))))
            if np.isfinite(v.real) and abs(v.real) > 0 and R >= 1.0:
                realdef = max(realdef, abs(v.imag) / abs(v.real))
            if (-1) ** n * v.real < -1e-7 * scale:
                fails.append(f"n={n}@R={R:g}")
    return fails, realdef


def main():
    print("== (1) CM scan, W(R)=exp(-c R^alpha), c=1, n<=8 ==")
    print("   R window [0.1, min(100,17^(1/alpha))]: below 0.1 the Cauchy formula for n=8 loses all digits,")
    print("   above it W<1e-8 and f^(n) is not resolvable in double precision")
    print("  alpha   R-range          sign-fails (-1)^n W^(n)<0       realness   verdict")
    for a in (0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 2.0, 3.0):
        Rhi = min(100.0, 17.0 ** (1 / a))
        Rg = np.logspace(-1, math.log10(Rhi), 21)
        f = (lambda a: (lambda z: np.exp(-np.asarray(z).astype(complex) ** a)))(a)
        fails, rd = cm_scan(f, Rg)
        tag = "IN (CM)" if not fails else "OUT"
        show = "-" if not fails else " ".join(fails[:2])
        print(f"  {a:4.1f}   [{Rg[0]:.1e},{Rg[-1]:.1e}]  {show:36s} {rd:8.1e}   {tag}")

    print("\n== (2) witness alpha=1/2: Levy(1/2) density, exact identity ==")
    u = np.linspace(-18, 18, 400001)
    t = np.exp(u)
    dens = (1 / (2 * np.sqrt(np.pi))) * t ** -1.5 * np.exp(-1 / (4 * t))
    mid = (t > 1e-2) & (t < 1e2)
    print(f"   density closed form = t^-3/2 exp(-1/4t)/(2 sqrt pi): min on t in [1e-2,1e2] = {dens[mid].min():.3e} > 0")
    print(f"   Int dens dt = {np.trapezoid(dens * t, u):.12f}")
    err = 0.0
    for R in (0.1, 0.5, 1.0, 4.0, 25.0):
        W = np.trapezoid(np.exp(-t * R) * dens * t, u)
        err = max(err, abs(W - np.exp(-np.sqrt(R))) / np.exp(-np.sqrt(R)))
    print(f"   max rel err  |Int e^-tR p(t)dt - e^-sqrt R|  over R in {{.1,.5,1,4,25}} = {err:.1e}")
    print("   => e^-sqrt(R) IS a Bernstein profile of a POSITIVE measure: 0<W<=1, no power law")

    print("\n== (3) engine: x^alpha = (alpha/Gamma(1-alpha)) Int_0^inf (1-e^-xt) t^-1-alpha dt ==")
    v = np.linspace(-200, 40, 2400001)
    s = np.exp(v)
    for a in (0.3, 0.5, 0.7):
        x = 1.7
        I = np.trapezoid((1 - np.exp(-x * s)) * s ** (-1 - a) * s, v)
        lhs = a / math.gamma(1 - a) * I
        print(f"   alpha={a}: formula={lhs:.10f}  x^alpha={x**a:.10f}  relerr={abs(lhs-x**a)/x**a:.1e}"
              "   (kernel (1-e^-xt)t^-1-a >= 0 => Bernstein function for alpha in (0,1))")

    print("\n== (4) alpha>1: W'' = c a R^(a-2) e^-cR^a (c a R^a - (a-1)), c=1 ==")
    for a in (1.2, 1.5, 2.0, 3.0):
        Rs = ((a - 1) / a) ** (1 / a)
        grid = 10.0 ** np.arange(-3, 1.01, 0.05)
        W2 = a * grid ** (a - 2) * np.exp(-grid ** a) * (a * grid ** a - (a - 1))
        print(f"   alpha={a}: W''<0 on (0,R*), R*=((a-1)/a)^(1/a)={Rs:.6f}; min W''={W2.min():+.3e} => NOT CM")

    print("\n== (5) log-power, EXACT profile: W = exp(-c (Log z)^beta) ==")
    print("   (a) branch: (Log z)^beta is non-real across R in (0,1) for beta not an integer")
    for b in (0.3, 0.5, 0.7, 1.5, 2.5):
        R, eps = 0.5, 1e-6
        vp = np.exp(-(np.log(complex(R, eps)) ** b))
        vm = np.exp(-(np.log(complex(R, -eps)) ** b))
        print(f"   beta={b:4.1f}: Im W(R+i0)={vp.imag:+.3e}  Im W(R-i0)={vm.imag:+.3e}"
              f"  Re-jump={abs(vp.real-vm.real):.3e}  -> real-valued CM W cannot match")
    print("   (b) integer beta: endpoint limit vs W(0+) = mu(R_+) in [0,1]")
    print("       (entries under e^-709 are clipped; the point is the limit, not the value)")
    for b in (1, 2, 3, 4):
        lo = math.exp(max(-709.0, min(709.0, -(math.log(1e-6)) ** b)))
        hi = math.exp(max(-709.0, min(709.0, -(math.log(1e6)) ** b)))
        why = ("unbounded as R->0 => |W|>1" if lo > 1 else
               ("W->0 as R->0 => mu=0 => W=0 (contradiction)" if lo < 1e-6 else "ok"))
        print(f"   beta={b}: W(1e-6)={lo:.3e}  W(1e6)={hi:.3e}   -> {why}")

    print("\n== (6) log-power, ASYMPTOTIC profile: rho(t)=t^-1 exp(-c|ln t|^beta), c=1 ==")
    print("   parametrise R=e^L; g(u)=-e^(L+u)-|u|^beta, ln W = max g + log Int exp(g-max g) du")
    for b in (0.5, 1.0, 2.0, 3.0):
        rows, vals = [], []
        for L in (10.0, 30.0, 100.0, 300.0, 600.0):
            u = np.linspace(-L - 40, min(40.0, -L + 690), 400001)
            A = L + u
            g = -np.exp(np.minimum(A, 700.0)) - np.abs(u) ** b
            g = np.where(A < 690, g, -np.inf)
            M = g.max()
            lnW = M + math.log(np.trapezoid(np.exp(g - M), u))
            rows.append(-lnW / L ** b)
            vals.append(-lnW)
        beff = (math.log(vals[-1]) - math.log(vals[-3])) / math.log(600.0 / 100.0)
        print(f"   beta={b}: -lnW/L^beta at L=10,30,100,300,600: " + " ".join(f"{x:.4f}" for x in rows)
              + f"   -> c=1; local exponent beta_eff(L=100..600) = {beff:.3f}")

    print("\n== (7) R^-p (ln R)^q ==")
    print("   q>0: d/dR[R^-p L^q] = R^-p-1 L^q-1 (q - p L) changes sign at ln R=q/p")
    for p, q in ((2.0, 1.0), (2.0, 2.0)):
        print(f"     p={p},q={q}: sign change at ln R={q/p:.3f} => not monotone => NOT CM")
    print("   q<0 (necessary-condition scan n<=8; R>=2 keeps the branch point R=1 outside the circle):")
    for p, q in ((0.0, -1.0), (0.5, -1.0), (1.0, -2.0), (2.0, -0.5)):
        f = (lambda p, q: (lambda z: np.asarray(z).astype(complex) ** (-p)
                           * np.log(np.asarray(z).astype(complex)) ** q))(p, q)
        fails, rd = cm_scan(f, 10.0 ** np.arange(0.3, 3.01, 0.2), rel=0.05)
        print(f"     p={p:3.1f} q={q:+.1f}: sign-fails={'-' if not fails else ' '.join(fails[:2]):20s}"
              f" realness={rd:.1e}")

    print("\n== (8) allowed set: gapped and gapless members both present ==")
    print("   member                   W(0+)   max W   p_eff=-dlnW/dlnR@[10,100]   class")
    Rg2 = np.array([10.0, 100.0])
    members = [("delta_m (m=0.7)", lambda R: math.exp(-0.7 * R), "gapped  W<=e^-0.7R"),
               ("U[1,2] uniform", lambda R: (math.exp(-R) - math.exp(-2 * R)) / R, "gapped  W<=e^-R"),
               ("stable alpha=0.5", lambda R: math.exp(-math.sqrt(R)), "no gap  exp(-sqrt R)"),
               ("stable alpha=0.3", lambda R: math.exp(-R ** 0.3), "no gap  exp(-R^0.3)"),
               ("stable alpha=0.9", lambda R: math.exp(-R ** 0.9), "no gap  stretched exp")]
    for name, W, cls in members:
        w0 = W(1e-8)
        mx = max(W(10.0 ** x) for x in np.arange(-8, 4, 0.1))
        pe = -(math.log(W(Rg2[1])) - math.log(W(Rg2[0]))) / math.log(Rg2[1] / Rg2[0])
        print(f"   {name:22s} {w0:7.4f}  {mx:6.4f}   {pe:20.4f}   {cls}")
    print("   => both a gapped and a gapless member: the two constraints have zero polarity on the gap")

    print("\n== (9) gap criterion: rate = inf supp rho ==")
    for a, b in ((0.3, 1.0), (0.7, 2.0)):
        R = 10.0 ** np.arange(0, 3.01, 1.0)
        W = (np.exp(-a * R) - np.exp(-b * R)) / R
        print(f"   rho=U[{a},{b}]:  R*W(R)/e^-{a}R  = " + " ".join(f"{x:.4f}" for x in (R * W / np.exp(-a * R)))
              + "   -> 1: rate = a = inf supp rho")
    print("   rho=U[1e-4,1]: R*W(R) at R=1e4 = "
          f"{1e4 * (math.exp(-1) - math.exp(-1e4)) / 1e4:.3e} vs e^-1 = {math.exp(-1):.3e}"
          "  -> rate 1e-4: mass near 0 sets the rate")


main()
