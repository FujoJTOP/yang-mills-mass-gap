"""kt_c_geom.py -- KT-RHO2 T1/T2: geometry, multiplicity patterns, double well.

C1  metric graphs (Kirchhoff vertices, Dirichlet boundary vertices) -- exact continuum
    Laplacians, discretised with a UNIFORM mesh width so vertex coupling is the true
    Kirchhoff one.  Intervals (DD/ND/NN), stars, random trees, lasso, figure eight.
C2  genuine 2D Dirichlet domains (5-point sparse FD, ghost-zero boundary): square,
    rectangle, L, dumbbell, annulus, straight tube, Y-tube.
C3  the double well: rho^2 = mu2/mu1 with the VACUUM SUBTRACTED (mu0 = 0), plus the
    bottom near-degeneracy (Kramers/tunnelling pattern).
Single process, deterministic.  Writes kt_c_geom.json + kt_c_geom.log per block.
"""
import json, numpy as np
from scipy.sparse import coo_matrix, identity
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh_tridiagonal

JSON, LOG = "kt_c_geom.json", []
OUT = {}


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.append(s)


def flush():
    with open(JSON, "w", encoding="utf-8") as f:
        json.dump(OUT, f, indent=1)
    with open("kt_c_geom.log", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


# ------------------------------------------------------------------ metric graphs
def graph_spec(edges, dirichlet, h=0.005, k=6):
    """-d^2/dx^2 on a metric graph, Kirchhoff vertices, Dirichlet nodes eliminated."""
    nid = {}
    key = []

    def newnode(tag):
        nid[tag] = len(key)
        key.append(tag)
        return nid[tag]

    # vertex tags are the given u,v; interior tags are (u,v,i)
    nid = {}
    key = []
    rows, cols, vals = [], [], []
    for (u, v, L) in edges:
        if u not in nid:
            newnode(u)
        if v not in nid:
            newnode(v)
        a, b = nid[u], nid[v]
        n = max(1, int(round(L / h)))
        hh = L / n
        prev = a
        for i in range(1, n):
            w = newnode((u, v, i))
            c = 1.0 / hh ** 2
            rows += [prev, w, prev, w]
            cols += [prev, w, w, prev]
            vals += [c, c, -c, -c]
            prev = w
        c = 1.0 / hh ** 2
        rows += [prev, b, prev, b]
        cols += [prev, b, b, prev]
        vals += [c, c, -c, -c]
    N = len(key)
    A = coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsc()
    keep = [i for i in range(N) if key[i] not in dirichlet]
    A = A[keep][:, keep].tocsc()
    lam = np.sort(eigsh(A, k=min(k, A.shape[0] - 2), sigma=-1e-4, which="LM",
                          return_eigenvectors=False))
    return lam


def group(v, tol=1e-2):
    vals, mult = [], []
    for z in np.sort(v):
        if vals and abs(z - vals[-1]) <= tol * max(1e-12, abs(vals[-1])):
            mult[-1] += 1
        else:
            vals.append(z); mult.append(1)
    return np.array(vals), np.array(mult)


P("=" * 100)
P("C1  metric graphs (continuum Laplacian, Kirchhoff vertices, h=0.005)")
P("=" * 100)
g = {}
cases = {
    "interval_DD": ([(0, 1, 1.0)], {0, 1}),
    "interval_ND": ([(0, 1, 1.0)], {1}),
    "interval_NN": ([(0, 1, 1.0)], set()),
    "interval_DD_len2": ([(0, 1, 2.0)], {0, 1}),
}
for m in [2, 3, 4, 5]:
    cases[f"star{m}_equal"] = ([(0, i + 1, 1.0) for i in range(m)], {i + 1 for i in range(m)})
cases["star3_unequal"] = ([(0, 1, 1.0), (0, 2, 1.0), (0, 3, 2.0)], {1, 2, 3})
cases["star3_oneD"] = ([(0, 1, 1.0), (0, 2, 1.0), (0, 3, 2.0)], {1, 2})
cases["lasso_loop2_tail1"] = ([(0, 1, 2 / 3), (1, 2, 2 / 3), (2, 0, 2 / 3), (0, 3, 1.0)], {3})
cases["lasso_loop6_tail1"] = ([(0, 1, 2.0), (1, 2, 2.0), (2, 0, 2.0), (0, 3, 1.0)], {3})
cases["fig8"] = ([(0, 1, 1.0), (1, 0, 1.0), (0, 2, 1.0), (2, 0, 1.0)], set())
for name, (edges, dr) in cases.items():
    lam = graph_spec(edges, dr, h=0.004, k=7)
    if abs(lam[0]) < 1e-6:          # constant mode = vacuum
        lam = lam[1:]
    vals, mult = group(lam, tol=1e-2)
    if len(vals) < 2:
        P(f"  {name:<20s} too few levels")
        continue
    d = dict(raw=float(lam[1] / lam[0]), merged=float(vals[1] / vals[0]),
             mult=[int(m) for m in mult[:5]], split=float(abs(lam[1] - lam[0]) / lam[0]),
             lam=[float(z) for z in lam[:5]])
    g[name] = d
    P(f"  {name:<20s} lam={np.round(lam[:5], 5)}  raw={d['raw']:.6f} merged={d['merged']:.6f}"
      f"  mult={d['mult']}")
OUT["C1_graphs"] = g
flush()

P("  random trees (3-5 arms, lengths U[0.3,1.6]), 300 samples, h=0.01")
rng = np.random.default_rng(7)
rs, worst = [], None
for _ in range(300):
    m = int(rng.integers(3, 6))
    edges = [(0, i + 1, float(rng.uniform(0.3, 1.6))) for i in range(m)]
    lam = graph_spec(edges, {i + 1 for i in range(m)}, h=0.01, k=5)
    vals, mult = group(lam, tol=1e-2)
    if len(vals) < 2:
        continue
    r = vals[1] / vals[0]
    rs.append(r)
    if worst is None or r > worst[0]:
        worst = (r, m, [round(e[2], 3) for e in edges])
rs = np.array(rs)
P(f"    sup merged rho^2 = {rs.max():.6f}   median = {np.median(rs):.6f}")
P(f"    worst: rho^2={worst[0]:.6f} m={worst[1]} lengths={worst[2]}")
OUT["C1_random_trees"] = dict(sup=float(rs.max()), median=float(np.median(rs)),
                              worst=dict(rho2=float(worst[0]), m=int(worst[1]), lens=worst[2]))
flush()


# ------------------------------------------------------------------ 2D domains
def domain_rho2(mask_fn, n=140, box=((0, 1), (0, 1)), k=5):
    (x0, x1), (y0, y1) = box
    xs = np.linspace(x0, x1, n)
    ys = np.linspace(y0, y1, n)
    hx, hy = xs[1] - xs[0], ys[1] - ys[0]
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    M = mask_fn(X, Y)
    idx = -np.ones(M.shape, int)
    idx[M] = np.arange(int(M.sum()))
    N = int(M.sum())
    rows, cols, vals = [], [], []
    for (a, b, coef) in ((1, 0, -1 / hx ** 2), (-1, 0, -1 / hx ** 2),
                         (0, 1, -1 / hy ** 2), (0, -1, -1 / hy ** 2)):
        js = list(range(max(0, -a), M.shape[0] - max(0, a)))
        ks = list(range(max(0, -b), M.shape[1] - max(0, b)))
        src = idx[np.ix_(js, ks)]
        dst = idx[np.ix_([j + a for j in js], [k2 + b for k2 in ks])]
        ok = (src >= 0) & (dst >= 0)
        rows += list(src[ok]); cols += list(dst[ok]); vals += [coef] * int(ok.sum())
    rows += list(np.arange(N)); cols += list(np.arange(N))
    vals += list((2 / hx ** 2 + 2 / hy ** 2) * np.ones(N))
    A = coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsc()
    lam = np.sort(eigsh(A, k=k, sigma=-1e-2, which="LM", return_eigenvectors=False))
    return lam


P("=" * 100)
P("C2  2D Dirichlet domains (5-point FD, ghost-zero boundary)")
P("=" * 100)
doms = {
    "square": (lambda X, Y: np.ones_like(X, bool), ((0, 1), (0, 1)), 120),
    "rect2to1": (lambda X, Y: np.ones_like(X, bool), ((0, 2), (0, 1)), 160),
    "Lshape": (lambda X, Y: ~((X > 0.5) & (Y > 0.5)), ((0, 1), (0, 1)), 140),
    "dumbbell": (lambda X, Y: (((X < 0.35) & (Y > 0.3) & (Y < 0.7)) |
                               ((X > 0.65) & (Y > 0.3) & (Y < 0.7)) |
                               ((X >= 0.35) & (X <= 0.65) & (Y > 0.47) & (Y < 0.53))),
                 ((0, 1), (0, 1)), 220),
    "annulus": (lambda X, Y: (((X - 0.5) ** 2 + (Y - 0.5) ** 2 < 0.25 ** 2) &
                              ((X - 0.5) ** 2 + (Y - 0.5) ** 2 > 0.12 ** 2)),
                ((0, 1), (0, 1)), 240),
    "tube_straight": (lambda X, Y: (Y > 0.45) & (Y < 0.55), ((0, 2), (0, 1)), 160),
    "Y_tube": (lambda X, Y: (((np.abs(X - 0.5) < 0.03) & (Y > 0.25) & (Y < 0.95)) |
                             ((np.abs(Y - 0.25) < 0.03) & (X > 0.1) & (X < 0.9))),
               ((0, 1), (0, 1)), 260),
}
d2 = {}
for name, (mf, box, n) in doms.items():
    lam = domain_rho2(mf, n=n, box=box, k=5)
    d2[name] = dict(lam=[float(z) for z in lam[:4]], rho2=float(lam[1] / lam[0]))
    P(f"  {name:<14s} n={n}  lam={np.round(lam[:4], 4)}  rho^2 = {lam[1]/lam[0]:.6f}")
P("  analytic: square 2.5 ; rect 2:1 = 3.4 ; thin tube -> 4")
OUT["C2_domains"] = d2
flush()

# ------------------------------------------------------------------ double well
P("=" * 100)
P("C3  double well V=Lam(x^2-1)^2 on [-3,3]: mu_n = E_n - E_0 (vacuum subtracted)")
P("=" * 100)
N = 8000
x = np.linspace(-3, 3, N)
h = x[1] - x[0]
E = np.full(N, 2.0 / h ** 2)
e = np.full(N - 1, -1.0 / h ** 2)
dw = {}
prev_rho = None
for Lam in [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 50.0]:
    V = Lam * (x ** 2 - 1) ** 2
    ev = eigh_tridiagonal(E + V, e, eigvals_only=True, select="i", select_range=(0, 6))
    mu = np.diff(ev[:6], prepend=ev[0])[1:]        # mu_n = E_n - E_0 for n>=1
    mu[0] = ev[1] - ev[0]
    mu = ev[1:] - ev[0]
    rho2 = mu[1] / mu[0]
    d = dict(E=[float(z) for z in ev[:5]], mu=[float(z) for z in mu[:4]],
             rho2=float(rho2), pair_split=float((ev[1] - ev[0]) / (ev[2] - ev[0])))
    dw[f"Lam{Lam:g}"] = d
    P(f"  Lam={Lam:<5g} mu1=E1-E0={mu[0]:.3e}  mu2={mu[1]:.5f}  rho2={rho2:.4f}"
      f"   bottom pair: (E1-E0)/(E2-E0)={d['pair_split']:.3e}"
      f"   growth vs prev: {(rho2/prev_rho if prev_rho else 1):.2f}x")
    prev_rho = rho2
for eps in [0.0, 1e-2, 5e-2, 0.2]:
    Lam = 20.0
    V = Lam * (x ** 2 - 1) ** 2 + eps * x
    ev = eigh_tridiagonal(E + V, e, eigvals_only=True, select="i", select_range=(0, 5))
    mu = ev[1:] - ev[0]
    dw[f"tilt_eps{eps:g}"] = dict(mu=[float(z) for z in mu[:3]], rho2=float(mu[1] / mu[0]))
    P(f"  Lam=20 tilt={eps:<6g} mu1={mu[0]:.4e} mu2={mu[1]:.5f} rho2={mu[1]/mu[0]:.4f}")
OUT["C3_double_well"] = dw
flush()

P("=" * 100)
P("C4  tally")
P("=" * 100)
P(f"  planar Dirichlet domains : max rho^2 = {max(v['rho2'] for v in d2.values()):.6f}")
P(f"  metric graphs            : max merged rho^2 = {max(v['merged'] for v in g.values()):.6f}")
P(f"  random trees             : sup = {rs.max():.6f}")
P(f"  double well              : rho^2 {dw['Lam0.5']['rho2']:.3f} (Lam=0.5) ->"
  f" {dw['Lam50']['rho2']:.3f} (Lam=50)")
flush()
P("done. json ->", JSON)
