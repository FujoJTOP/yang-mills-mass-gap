"""S08-PSI: quantify psi = sum_{i<j} log d_ij^2  (the r12 correction to r08's Phi).

Question: is  Phi_O09 = sum g(d)  =  Phi_r08 + psi  with  Phi_r08 = -log|Delta|^2 ?
And what makes  E[Phi'^2] = E[Phi'']  collapse?

Model = r08/r12 Weyl reduction on the SU(3) torus:
  theta in T^3,  lam in t (sum lam_i = 0),  tau: theta -> theta + tau*lam
  |Delta|^2 = prod_{i<j} 4 sin^2(d_ij/2),  d_ij -> d_ij + tau*(lam_i-lam_j)   (AFFINE in tau)
  weight    = |Delta|^2 dtheta   (= e^{-Phi_r08}, the Gibbs measure of Phi_r08 itself)
Conformal coords (r12): d12=a, d23=b, d13=a+b,  a,b in (0,2pi).

numpy only (no sympy needed at runtime beyond the algebra block).
"""
import numpy as np
import sympy as sp

say_lines = []


def say(s=""):
    print(s)
    say_lines.append(s)


say("S08-PSI  |  psi = sum_{i<j} log d_ij^2   (role of the r12 correction)")
say("=" * 78)

# ---------------------------------------------------------------- 1. algebra
d, a, v = sp.symbols("d a v", positive=True)
g = -2 * sp.log(sp.sin(d / 2) / (d / 2))          # O09/q11 exponential-coord Jacobian
phi08 = -sp.log(4 * sp.sin(d / 2) ** 2)           # r08:  -log|Delta|^2  (per pair)
psi1 = sp.log(d ** 2)                             # psi (per pair)

diff = sp.logcombine(g - phi08 - psi1, force=True)
say("[1] per-pair  g(d) - (-log 4sin^2(d/2)) - log d^2  =  %s   (simplify: %s)"
    % (diff, sp.simplify(diff)))
say("    numeric d=1.234 : %.3e   (=> Phi_O09 - Phi_r08 = sum log d^2 exactly, no constant)"
    % float((g - phi08 - psi1).subs(d, 1.2343).evalf()))

gp, gpp = sp.simplify(sp.diff(g, d)), sp.simplify(sp.diff(g, d, 2))
say("[2] g'  = %s        [= -cot(d/2) + 2/d :  %s]"
    % (gp, sp.simplify(gp + sp.cot(d / 2) - 2 / d)))
say("    g'' = %s   [= csc^2(d/2)/2 - 2/d^2 : %s]"
    % (gpp, sp.simplify(gpp - sp.csc(d / 2) ** 2 / 2 + 2 / d ** 2)))

D = d + a * v
psiv = sp.log(D ** 2)
say("[3] psi(v) = log[(d+a v)^2] :  psi' = %s ,  psi'' = %s"
    % (sp.simplify(sp.diff(psiv, v)), sp.simplify(sp.diff(psiv, v, 2))))
f1, f2, f3 = (sp.simplify(sp.diff(phi08, d, k)) for k in (1, 2, 3))
say("[4] f(d) = -log4sin^2(d/2):  f'=%s  f''=%s  f'''=%s   (f''' != 0 => Phi_r08 NOT quadratic in v)"
    % (f1, f2, f3))

f3f = sp.lambdify(d, f3, "math")
d3 = sum(av ** 3 * f3f(dv) for av, dv in [(1.0, 1.1), (1.0, 1.7), (2.0, 2.8)])
say("    Phi_r08''' at generic (d;a) = (1.1;1),(1.7;1),(2.8;2):  %.6f  (nonzero)" % d3)

# ------------------------------------------------- 2. quadrature on T^2 (r08)
def ensemble(N):
    dx = 2 * np.pi / N
    x = (np.arange(N) + 0.5) * dx
    A = x[:, None]
    B = x[None, :]
    d12, d23, d13 = A, B, A + B
    w = (4 * np.sin(d12 / 2) ** 2) * (4 * np.sin(d23 / 2) ** 2) * (4 * np.sin(d13 / 2) ** 2)
    return w, [d12, d23, d13], dx


def averages(lam, N):
    w, ds, dx = ensemble(N)
    Z = w.sum() * dx * dx
    with np.errstate(divide="ignore", invalid="ignore"):
        av = np.array([lam[0] - lam[1], lam[1] - lam[2], lam[0] - lam[2]], float)
        V = sum(-np.log(4 * np.sin(dp / 2) ** 2) for dp in ds)            # Phi_r08
        Vp = sum(-av[i] / np.tan(ds[i] / 2) for i in range(3))            # Phi_r08'
        Vpp = sum(0.5 * av[i] ** 2 / np.sin(ds[i] / 2) ** 2 for i in range(3))
        ps = sum(np.log(dp ** 2) for dp in ds)                            # psi
        psp = sum(2 * av[i] / ds[i] for i in range(3))                    # psi'
        pspp = sum(-2 * av[i] ** 2 / ds[i] ** 2 for i in range(3))        # psi''
    E = lambda F: float((w * F).sum() * dx * dx / Z)
    out = dict(Z=Z, Vpp=E(Vpp), Vp2=E(Vp * Vp), Vp=E(Vp), pspp=E(pspp),
               psp2=E(psp * psp), psp=E(psp), Phi=E(V + ps), Phipp=E(Vpp + pspp),
               Phip2=E((Vp + psp) ** 2), Vpsp=E(Vp * psp), Phippsp=E((Vp + psp) * psp),
               a2d2=E(sum(av[i] ** 2 / ds[i] ** 2 for i in range(3))), N=N, lam=lam)
    return out


say("")
say("[5] midpoint quadrature, weight |Delta|^2 dtheta  (r08 measure = Gibbs measure of Phi_r08)")
say("    %-12s %4s %10s %10s %10s %10s %10s" % ("lam", "N", "E[V'']", "E[V'^2]", "S-D(O09)", "B=-Eps''-Varps'", "C=-Cov"))
for lam in [(1, 0, -1), (1, -1, 0)]:
    for N in (1024, 1536):
        o = averages(lam, N)
        EPhi_p = o["psp"]                                    # E[V'] = 0  =>  E[Phi'] = E[psi']
        varPhi_p = o["Phip2"] - EPhi_p ** 2                  # Var(Phi')
        SD = o["Phipp"] - varPhi_p
        B = -o["pspp"] - (o["psp2"] - o["psp"] ** 2)
        C = -(o["Phippsp"] - EPhi_p * o["psp"])
        say("    %-12s %4d %10.6f %10.6f %10.6f %10.6f %10.6f"
            % (str(lam), N, o["Vpp"], o["Vp2"], SD, B, C))
say("")
for lam in [(1, 0, -1), (1, -1, 0)]:
    o = averages(lam, 1536)
    EPhi_p = o["psp"]
    varpsip = o["psp2"] - o["psp"] ** 2
    varPhi_p = o["Phip2"] - EPhi_p ** 2
    SD = o["Phipp"] - varPhi_p
    B = -o["pspp"] - varpsip
    C = -(o["Phippsp"] - EPhi_p * o["psp"])
    say("    lam=%s : E[psi']=%+.6f  E[psi'']=%+.6f  Var(psi')=%.6f  sum E[a^2/d^2]=%.6f (= -E[psi'']/2)"
        % (str(lam), o["psp"], o["pspp"], varpsip, o["a2d2"]))
    say("              E[Phi'']=%.6f  E[Phi'^2]=%.6f  Var(Phi')=%.6f  E[V']=%.2e"
        % (o["Phipp"], o["Phip2"], varPhi_p, o["Vp"]))
    say("              id.resid E[V'^2]-E[V'']=%+.3e | broken by E[Phi'^2]-E[Phi'']=%+.6f"
        % (o["Vp2"] - o["Vpp"], o["Phip2"] - o["Phipp"]))
    say("              S-D: direct=%+.6f  -E[psi'']-Var(psi')=%+.6f  -Cov(Phi',psi')=%+.6f"
        % (SD, B, C))
    say("")

# exact-rational check of the CT statement (r08 sec 3.1) is already in r08; here we
# confirm the two ATTRIBUTIONS of the identity, not its value.
say("[6] attribution test:  weight w = |Delta|^2 = e^{-Phi_r08}  =>  int d_tau^2 w dtheta = 0")
say("    (translation invariance)  <=>  E[Phi_r08'^2] = E[Phi_r08'']  for ANY smooth Phi_r08,")
say("    whether or not it is quadratic in tau.  Numeric: see column 3 vs 4 above (they agree).")
say("")
print("done")
