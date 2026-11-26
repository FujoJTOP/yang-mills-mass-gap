"""
O09 - Phi'' provenance and magnitude.  py -3.12, numpy only, deterministic quadrature.

Object (docs/paper-v2/ch05-theorems.tex:196-240, ch16-sixth-species.tex:120-135):
    flat chart X = log U in su(N),  dU = J dX,  Phi = -log J
    SU(N):  Phi = sum_{i<j} g(mu_i - mu_j),  g(d) = -2 log( sin(d/2)/(d/2) )
            equivalently Phi = -log |c_SU(N)/c_flat|^2 ,  c_flat = prod_{a>0} alpha(X)
    Phi'' := <v, Hess_X Phi v>,  Phi' := <v, grad_X Phi>,  xi = <v,X>
    G''(a) = E[H_vv + Phi''|xi=a] - Var(G_v + Phi'|xi=a),   H := Hess S_W
      (eq:split adds Phi'' separately; using H := Hess V would double count Phi'')
    normalization tr(T^a T^b)=delta_ab/2, v = T_3 = diag(1/2,-1/2,0..), <v,v> := 2 tr v^2 = 1
    -> Hess Phi|_0 = (N/12) <v,v> = N/12  (1/6, 1/4, 1/3, ... for N=2,3,4)
"""
import sys
import numpy as np

LOG = open(r"D:\Dev\math\w7\o09_phicurv.log", "w", encoding="utf-8")


def say(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True)
    LOG.write(line + "\n")
    LOG.flush()


np.set_printoptions(precision=10, suppress=True)
TWO_PI = 2 * np.pi


# ---------------- g and derivatives (series-safe at d -> 0) ----------------
def g(d):
    d = np.asarray(d, float)
    small = np.abs(d) < 1e-4
    safe = np.where(small, 1.0, d)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = -2.0 * np.log(np.sin(safe / 2) / (safe / 2))
    ser = d ** 2 / 12 + d ** 4 / 720 + d ** 6 / 362880
    return np.where(small, ser, out)


def phi_from_mu(mu):
    n = mu.shape[-1]
    i, j = np.triu_indices(n, 1)
    return g(mu[..., i] - mu[..., j]).sum(axis=-1)


# ---------------- cached Gauss-Legendre on [0,1] ----------------
GL = {}


def gl(n):
    if n not in GL:
        x, w = np.polynomial.legendre.leggauss(n)
        GL[n] = (0.5 * (x + 1.0), 0.5 * w)      # nodes/weights on [0,1]
    return GL[n]


say("=" * 100)
say("[1] g''(d) = csc^2(d/2)/2 - 2/d^2 = 1/6 + d^2/120 + ...  on |d| < 2pi")
dd = np.linspace(1e-3, TWO_PI - 1e-2, 400001)
gg2 = 0.5 / np.sin(dd / 2) ** 2 - 2.0 / dd ** 2
say(f"  series check  g''(0) = 1/6 exactly (code series 1/6 + d^2/120)")
say(f"  min g'' on (0, 2pi-1e-2] = {gg2.min():.12f} at d = {dd[gg2.argmin()]:.6f}")
m1 = dd <= np.pi
m2 = dd > np.pi
say(f"  non-increasing steps: on (0,pi] = {int((np.diff(gg2[m1]) < 0).sum())}, "
    f"on (pi,2pi) = {int((np.diff(gg2[m2]) < 0).sum())}  -> g'' strictly increasing on each half")
say(f"  floor: g'' >= 1/6 = {1/6:.12f} on |d| < 2pi, equality only as d -> 0")
say(f"  g''(2pi-0.01) = {gg2[-1]:.6f}  (diverges as d->2pi^-, but J ~ (2pi-d)^2 kills the integrand)")
say(f"  => Phi'' = sum_ij g''(mu_ij) (dmu_ij)^2 >= (1/6) sum_ij (dmu_ij)^2  pointwise, any X")

say("=" * 100)
say("[2] Hess Phi|_0 along v = T_3 : expect N/12  (floors 1/6, 1/4, 1/3 ...)")
for N in (2, 3, 4, 5, 6):
    mu = np.zeros(N)
    mu[0], mu[1] = 0.5, -0.5                    # eigenvalues of v = T_3
    h = 1e-5
    F = lambda t: float(phi_from_mu(t * mu))
    fpp = (-F(2 * h) + 16 * F(h) - 30 * F(0) + 16 * F(-h) - F(-2 * h)) / (12 * h ** 2)
    say(f"  N={N}: Hess Phi|_0 = {fpp:.10f}   N/12 = {N/12:.10f}   diff = {fpp-N/12:+.2e}")
say("  (ch16-sixth-species.tex:130-135 measured 1/6, 1/4, 1/3 for N=2,3,4 -> same numbers)")

say("=" * 100)
say("[3] pointwise Phi'' at X != 0 : random Hermitian X in su(N), random traceless unit v")


def su_gens(N):
    G = []
    for i in range(N):
        for j in range(i + 1, N):
            A = np.zeros((N, N), complex); A[i, j] = 1; A[j, i] = 1; G.append(A / 2)
            B = np.zeros((N, N), complex); B[i, j] = -1j; B[j, i] = 1j; G.append(B / 2)
    for k in range(1, N):
        D = np.zeros((N, N), complex)
        for i in range(k):
            D[i, i] = 1.0
        D[k, k] = -k
        D /= np.sqrt(2 * k * (k + 1))
        G.append(D)
    return G


rng = np.random.default_rng(7)
for N in (2, 3, 4):
    G = su_gens(N)
    vals = []
    for _ in range(30):
        X = sum([rng.normal() * T for T in G], np.zeros((N, N), complex))
        v = sum([rng.normal() * T for T in G], np.zeros((N, N), complex))
        v = v / np.sqrt(2 * np.real(np.trace(v @ v)))
        h = 1e-4
        F = lambda e: float(phi_from_mu(np.linalg.eigvalsh(X + e * v)))
        vals.append((-F(2 * h) + 16 * F(h) - 30 * F(0) + 16 * F(-h) - F(-2 * h)) / (12 * h ** 2))
    vals = np.array(vals)
    say(f"  N={N}: min Phi'' over 30 random (X,v) = {vals.min():.6f}  max = {vals.max():.4f}  "
        f"floor N/12 = {N/12:.6f}  -> {'all >= floor' if vals.min() >= N/12 - 1e-6 else 'SOME BELOW FLOOR'}")

# =====================================================================================
# [4] SU(2) single plaquette: cylindrical quadrature on the slice <v,xi> = a
#     xi in R^3 = su(2); v = e_3; rho=|xi|; u = xi_3/rho; J = sin^2(rho/2)/rho^2
#     S = beta (1 - cos(rho/2)cos(th/2) + sin(rho/2)sin(th/2) u),  Phi = -log J
#     d^3xi = 2pi s ds dxi_3  ->  slice weight  e^{-(S+Phi)} s ds
# =====================================================================================
say("=" * 100)
say("[4] SU(2) single plaquette, beta = 2.5 : Phi''(a) curve and the H1 ledger")
RMAX = TWO_PI - 1e-3
NS = 2000
HH = 1e-4


def slice_a(a, beta, theta):
    if abs(a) >= RMAX - 1e-5:
        return None
    smax = np.sqrt(max(RMAX ** 2 - a ** 2, 0.0))
    nd, wt = gl(NS)
    s = smax * nd
    wgt = smax * wt
    c, sth = np.cos(theta / 2), np.sin(theta / 2)

    def fields(x3):
        rho = np.sqrt(s ** 2 + x3 ** 2)
        rho = np.maximum(rho, 1e-12)
        u = x3 / rho
        sig = 1.0 - np.cos(rho / 2) * c + np.sin(rho / 2) * sth * u
        Phi = -2.0 * np.log(np.sin(rho / 2) / (rho / 2))
        return sig, Phi

    sig, Phi = fields(a)
    sp, Pp = fields(a + HH)
    sm, Pm = fields(a - HH)
    s2p, P2p = fields(a + 2 * HH)
    s2m, P2m = fields(a - 2 * HH)
    dsig = (sp - sm) / (2 * HH)
    dPhi = (Pp - Pm) / (2 * HH)
    ddsig = (-s2p + 16 * sp - 30 * sig + 16 * sm - s2m) / (12 * HH ** 2)
    ddPhi = (-P2p + 16 * Pp - 30 * Phi + 16 * Pm - P2m) / (12 * HH ** 2)
    w = np.exp(-(beta * sig + Phi)) * s * wgt
    Z = w.sum()
    m = lambda q: float((w * q).sum() / Z)
    Gv = beta * dsig
    Fv = Gv + dPhi
    out = dict(Phi2=m(ddPhi), H=m(beta * ddsig),
               varG=m(Gv ** 2) - m(Gv) ** 2,
               varPhi=m(dPhi ** 2) - m(dPhi) ** 2,
               cov=m(Gv * dPhi) - m(Gv) * m(dPhi),
               varF=m(Fv ** 2) - m(Fv) ** 2)
    out["Gpp"] = m(beta * ddsig + ddPhi) - out["varF"]
    return out


for th, tag, ref in ((np.pi, "theta=pi", (0.2051153, 0.0, 0.069121, 0.135995)),
                     (0.0, "theta=0 ", (0.176864, 0.553061, 0.0, 0.729924))):
    r = slice_a(0.0, 2.5, th)
    say(f"  VALIDATION a=0, beta=2.5, {tag} | this run vs H1A-EXACT.md §2.1:")
    say(f"    E[Phi'']  = {r['Phi2']:.7f}   ref {ref[0]:.7f}")
    say(f"    E[H_vv]   = {r['H']:+.3e}   ref {ref[1]:+.6f}")
    say(f"    Var(G)    = {r['varG']:.7f}   ref {ref[2]:.7f}")
    say(f"    G''(a=0)  = {r['Gpp']:+.7f}   ref {ref[3]:+.7f}")
    say(f"    Var(Phi') = {r['varPhi']:.3e}   Cov(G,Phi') = {r['cov']:+.3e}")

say("-" * 100)
say("  a-curve at beta=2.5 (lattice flat-chart units).  E[H_vv] = action Hessian only")
say("  %7s | %9s %9s | %8s %9s %8s | %9s | %10s" %
    ("a", "E[Phi2]", "E[H_vv]", "Var(G)", "Var(Phi1)", "Cov", "Gpp", "Phi2/VarG"))
a_grid = np.concatenate([np.linspace(0, 2.0, 9), np.linspace(2.25, 5.0, 6),
                         np.linspace(5.25, TWO_PI - 0.15, 4)])
curves = {}
for th in (np.pi, 0.0):
    say(f"  --- theta = {th/np.pi:.0f} pi ---")
    rows = []
    for a in a_grid:
        r = slice_a(float(a), 2.5, th)
        if r is None:
            continue
        rows.append((float(a), r))
        ratio = r["Phi2"] / r["varG"] if r["varG"] > 1e-13 else float("inf")
        say("  %7.3f | %9.5f %+9.5f | %8.5f %9.5f %+8.5f | %+9.5f | %10.2f" %
            (a, r["Phi2"], r["H"], r["varG"], r["varPhi"], r["cov"], r["Gpp"], ratio))
    curves[th] = rows

say("=" * 100)
say("[5] beta-scan: beta*(G''=0) and beta where Phi'' alone would fail (theta=pi)")
betas = np.array([0.0, 1.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.30661, 5.0, 6.0, 8.0])
for a in (0.0, 0.5, 1.0, 1.5, 2.0, 3.0):
    row = []
    for b in betas:
        r = slice_a(a, b, np.pi)
        row.append(f"b={b:.2f}:Pv={r['Phi2']:.4f},Vf={r['varF']:.4f},G''={r['Gpp']:+.4f}")
    say(f"  a={a:.2f}")
    for k in range(0, len(row), 3):
        say("      " + " | ".join(row[k:k + 3]))

    def root(fun, lo=0.0, hi=60.0, it=30):
        flo, fhi = fun(lo), fun(hi)
        if flo * fhi > 0:
            return None
        for _ in range(it):
            mid = 0.5 * (lo + hi)
            if fun(mid) * flo > 0:
                lo, flo = mid, fun(mid)
            else:
                hi = mid
        return 0.5 * (lo + hi)

    b1 = root(lambda b: slice_a(a, b, np.pi)["Gpp"])
    b2 = root(lambda b: slice_a(a, b, np.pi)["Phi2"] - slice_a(a, b, np.pi)["varF"])
    say(f"      beta*(G''=0) = {b1 if b1 is None else round(b1,5)}    "
        f"beta where Phi''=Var(force) = {b2 if b2 is None else round(b2,5)}")

say("=" * 100)
say("[6] does Phi'' alone cover Var(force) at beta=2.5, theta=pi ?")
worst = (1e9, None)
bad = []
for a, r in curves[np.pi]:
    q = r["Phi2"] - r["varF"]
    if q < 0:
        bad.append((round(a, 3), round(q, 5)))
    if q < worst[0]:
        worst = (q, a)
say(f"  min over a-grid of (Phi'' - Var(force)) = {worst[0]:+.5f} at a = {worst[1]:.3f}")
say(f"  a where Phi'' < Var(force): {bad if bad else 'none on the grid'}")
say("=" * 100)
LOG.close()
