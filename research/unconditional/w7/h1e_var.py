"""H1E — analytic Var(<v,grad S>|0) on one SU(2) plaquette; closed form of its beta^2 coefficient.

Question (YG): is delta(beta) = 0.2051153 - 0.0110593 beta^2 with the beta^2 coefficient
DERIVED (closed form) or merely fitted?  Model A convention = research/v3/species/w7/na_plaquette.py:
    U = exp(xi.sig/2), |xi|<2pi, h=rho/2, u=xi3/rho,
    S = beta(1 - cos h cos(th/2) + sin h sin(th/2) u),   Phi = -log J,  J = 4 sin^2 h / rho^2,
    v = e3, condition <v,xi> = a0 = 0  =>  xi3 = 0.
At th=pi: S|slice = beta (const) and G := <v,grad S> = beta sin(rho/2)/rho  exactly.
Slice measure: Lebesgue d^2xi times J  =>  density  sin^2(rho/2)/rho drho  (beta-independent).

Analytic integrals (x = rho/2 in (0,pi)):
    Z  = int_0^pi sin^2 x / x        dx = Cin(2pi)/2
    N1 = int_0^pi sin^3 x / x^2      dx = (3/4)[ln3 + Ci(pi) - Ci(3pi)]
    N2 = int_0^pi sin^4 x / x^3      dx = ln2 + Ci(2pi) - Ci(4pi)
         via IBP twice: d/dx[sin^4 x/x^2] = 4sin^3x cosx/x^2 - 2sin^4x/x^3,
         then d/dx[sin^3x cosx/x], integrand 3sin^2x cos^2x - sin^4x = (cos2x - cos4x)/2
    Var/beta^2 = N2/(4Z) - (N1/(2Z))^2 = (N2 Z - N1^2)/(4 Z^2)
Only Ci at pi, 2pi, 3pi, 4pi + logs appear.  No fitting.
"""
import numpy as np, math

EULER = 0.5772156649015328606
LN = math.log


def gl(f, a, b, n=2000):
    x, w = np.polynomial.legendre.leggauss(n)
    t = 0.5 * (b - a) * x + 0.5 * (a + b)
    return 0.5 * (b - a) * float(np.sum(w * f(t)))


def Ci_num(x):
    """Ci(x) = gamma + ln x + int_0^x (cos t - 1)/t dt  (integrand analytic, ~ -t/2 at 0)."""
    return EULER + LN(x) + gl(lambda t: (np.cos(t) - 1.0) / t, 0.0, x, 400)


def Cin_num(x):
    return EULER + LN(x) - Ci_num(x)


# ---------------------------------------------------------------- 1D analytic integrals
PI = math.pi
Z1 = gl(lambda x: np.sin(x) ** 2 / x, 0.0, PI, 4000)
N1 = gl(lambda x: np.sin(x) ** 3 / x ** 2, 0.0, PI, 4000)
N2 = gl(lambda x: np.sin(x) ** 4 / x ** 3, 0.0, PI, 4000)

# closed forms
Zc = 0.5 * Cin_num(2 * PI)
N1c = 0.75 * (LN(3.0) + Ci_num(PI) - Ci_num(3 * PI))
N2c = LN(2.0) + Ci_num(2 * PI) - Ci_num(4 * PI)

Ef = N1 / (2 * Z1)
Ef2 = N2 / (4 * Z1)
V0 = Ef2 - Ef * Ef
V0c = (N2c * Zc - N1c ** 2) / (4 * Zc ** 2)

print("=" * 92)
print("[1] slice integrals: quadrature vs closed form")
print(f"    Z  = {Z1:.15f}   closed 0.5*Cin(2pi) = {Zc:.15f}   d={Z1-Zc:+.3e}")
print(f"    N1 = {N1:.15f}   closed (3/4)[ln3+Ci(pi)-Ci(3pi)] = {N1c:.15f}   d={N1-N1c:+.3e}")
print(f"    N2 = {N2:.15f}   closed ln2+Ci(2pi)-Ci(4pi)   = {N2c:.15f}   d={N2-N2c:+.3e}")
print(f"    Ci(pi)={Ci_num(PI):.12f} Ci(2pi)={Ci_num(2*PI):.12f} Ci(3pi)={Ci_num(3*PI):.12f} Ci(4pi)={Ci_num(4*PI):.12f}")
print(f"    Cin(2pi) = {Cin_num(2*PI):.12f}")

print("\n[2] Var/beta^2  (theta=pi, slice = {xi3=0}, pure Haar)")
print(f"    E[G]/beta  = N1/(2Z) = {Ef:.12f}")
print(f"    E[G^2]/beta^2 = N2/(4Z) = {Ef2:.12f}")
print(f"    V0 (quadrature) = {V0:.15f}")
print(f"    V0 (closed form) = {V0c:.15f}   |d| = {abs(V0-V0c):.2e}")
ARCH = 0.0110593
print(f"    archived coefficient      = {ARCH:.7f}")
print(f"    |V0 - 0.0110593| / 0.0110593 = {abs(V0-ARCH)/ARCH:.3e}")

print("\n[3] Haar anchors of the same model (theta=pi)")
W0 = 1.0 / (2.0 * Cin_num(2 * PI))          # corpus anchor 1/(2Cin(2pi))
print(f"    1/(2Cin(2pi))            = {W0:.12f}   (corpus eq:v3anchor 0.20511530)")
print(f"    1/(4Z),  Z=Cin(2pi)/2     = {1.0/(4*Z1):.12f}")
print(f"    beta* = sqrt(W0/V0)       = {math.sqrt(W0/V0):.8f}   (corpus 4.30661)")
print(f"    V0*W0 = {V0*W0:.12f} ;  V0/W0 = {V0/W0:.12f} ;  1/(4 pi^2)*? = {1/(4*PI**2):.12f}")
for nm, cand in [("1/90", 1/90), ("1/(8pi^2)", 1/(8*PI**2)), ("(pi-3)/12", (PI-3)/12),
                 ("1/(2 pi^2 N^2),N=2", 1/(2*PI**2*4)), ("W0^2*1.07", W0*W0)]:
    print(f"    candidate {nm:20s} = {cand:.12f}   ratio V0/cand = {V0/cand:.6f}")


# ---------------------------------------------------------------- 2D slice quadrature (mimics na_plaquette.py)
def onb(v):
    v = np.asarray(v, float); v = v / np.linalg.norm(v)
    a = np.array([1.0, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(v, a); e1 /= np.linalg.norm(e1)
    return v, e1, np.cross(v, e1)


def moments(a0, v, beta, theta, ns=800, nph=64):
    v, e1, e2 = onb(v)
    R = math.sqrt(max(4 * PI ** 2 - a0 * a0, 1e-14))
    xs, ws = np.polynomial.legendre.leggauss(ns)
    xp, wp = np.polynomial.legendre.leggauss(nph)
    s = 0.5 * R * (xs + 1.0); wsv = 0.5 * R * ws
    ph = PI * (xp + 1.0); wpv = PI * wp
    S, PH = np.meshgrid(s, ph, indexing="ij")
    WS, WPH = np.meshgrid(wsv, wpv, indexing="ij")
    xi = (a0 * v + (S * np.cos(PH))[..., None] * e1 + (S * np.sin(PH))[..., None] * e2).reshape(-1, 3)
    w = (WS * S * WPH).reshape(-1)
    rho = np.linalg.norm(xi, axis=1)
    h = 0.5 * rho
    u = xi[:, 2] / rho
    Sv = beta * (1.0 - np.cos(h) * math.cos(theta / 2) + np.sin(h) * math.sin(theta / 2) * u)
    J = 4.0 * np.sin(h) ** 2 / rho ** 2                       # e^{-Phi}
    # G = <e3, grad S> with v = e3:  d/dxi3 of S, chain rule on rho and u
    dS_drho = beta * (0.5 * np.sin(h) * math.cos(theta / 2) + 0.5 * np.cos(h) * math.sin(theta / 2) * u)
    dS_du = beta * np.sin(h) * math.sin(theta / 2)
    # u = xi3/rho => d u/d xi3 = (1 - u^2)/rho ; d rho/d xi3 = u  (on the whole ball)
    G = dS_drho * u + dS_du * (1.0 - u * u) / rho
    lw = np.log(w) - Sv - (-np.log(J))
    mx = lw.max(); ww = np.exp(lw - mx); Zs = ww.sum()
    E = lambda f: float((ww * f).sum() / Zs)
    EG = E(G); return EG, E(G * G) - EG * EG, E(Sv)


print("\n[4] beta-strictness: Var(G)/beta^2 on the slice (2D quadrature, a0=0, v=e3)")
print(f"    {'theta':>8}{'beta':>7}{'E[S]':>10}{'Var(G)':>12}{'Var/beta^2':>14}")
for th in (PI, 0.9 * PI, 0.5 * PI):
    for beta in (0.5, 1.0, 2.5, 6.0):
        EG, V, ES = moments(0.0, (0, 0, 1.0), beta, th)
        print(f"    {th:8.4f}{beta:7.2f}{ES:10.4f}{V:12.6f}{V/beta**2:14.9f}")
    print()

print("[5] the same at theta=pi with a0 = 0 only: E[S] must be exactly beta (beta-independent slice)")
for beta in (0.5, 2.5, 6.0):
    EG, V, ES = moments(0.0, (0, 0, 1.0), beta, PI)
    print(f"    beta={beta:4.1f}  E[S]={ES:.12f}  beta={beta:.12f}  d={ES-beta:+.2e}   Var/beta^2={V/beta**2:.12f}")

print("\n[6] N-scaling input: operator-norm bound on the gradient of g=(1/N)Re tr(U B)")
print("    |d_v g| <= sqrt(tr(T_v T_v))/N * sqrt(tr(UU^dag)) = (1/sqrt2)/N * sqrt(N) = 1/sqrt(2N)")
for N in (2, 3, 4, 6, 10):
    print(f"    N={N:3d}  bound Var <= 1/(2N) = {1/(2*N):.6f}   |grad g| <= {1/math.sqrt(2*N):.6f}")
print("\nCLay 0%.  No QED.")
