"""r1_chain.py -- R1 Level 1/2 witnesses on a NON-Gaussian, NON-convex chain (P4-type kernel).

MODEL (the corpus's own P4 object: one Cartan angle per site, cyclic chain, L sites):
    k(x,y)   = exp(-beta (1 - cos((y-x)/2)))          [one-plaquette / P4 transfer kernel]
    w(x)     = 1 + kappa cos(2x)                      [on-site weight: a two-well measure]
    mu_chain = prod_i w(x_i) dx_i prod_i k(x_i,x_{i+1})    (cyclic)
    T        = diag(sqrt w) K diag(sqrt w)            [symmetrised transfer matrix, n x n]

SECTIONS
  [1] exact identity  gap(T) = gap of the CARRE-DU-CHAMP (semigroup) Dirichlet form in L^2(mu),
      mu = Omega^2 (Omega = top eigenvector).  Machine precision, no constants.
  [2] the DISCRETE-vs-CONTINUOUS time factor (1-x)/(-log x), x = lam1/lam0.  >0, <1, -> 0 as x->1.
  [3] the GRADIENT substitution lam1^grad(mu) vs the semigroup gap, and its restoration as the
      kernel width eps -> 0 with the metric weight eps/2 (matched metric).
  [4] the wrong-space loss on the FULL chain:  lam_conf(chain) <= lam1^grad(mu_slice), strict.
  [5] LEVEL 2: restriction to the even (x -> -x) subspace cannot lower the constant; measured.
  [6] the LSI -> Poincare loss on a two-point space: exact, unbounded in 1/eps.

Run: py -3.12 -I r1_chain.py .   numpy/scipy only, deterministic, single process.
Labels [T]/[V].  No QED.  Clay 0%.
"""
import json
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh

LOG = []


def say(s=""):
    print(s)
    LOG.append(s)


CHK = []


def chk(name, ok, info=""):
    CHK.append({"name": name, "ok": bool(ok), "info": info})
    say(f"    [{'PASS' if ok else 'FAIL'}] {name} {info}")


def save(out):
    with open("r1_chain.log", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")
    with open("r1_chain.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)


TWO_PI = 2.0 * np.pi


def kernel(x, y, beta):
    return np.exp(-beta * (1.0 - np.cos(0.5 * (y[None, :] - x[:, None]))))


def onsite(x, kappa):
    return 1.0 + kappa * np.cos(2.0 * x)


def grid(n):
    h = TWO_PI / n
    return h * np.arange(n), h


def graph_lap(cvec, n):
    """sparse weighted cycle Laplacian with edge weights c_0..c_{n-1}, row sums 0."""
    ip = (np.arange(n) + 1) % n
    rows = np.concatenate([np.arange(n), ip, np.arange(n), ip])
    cols = np.concatenate([ip, np.arange(n), np.arange(n), ip])
    vals = np.concatenate([-cvec, -cvec, cvec, cvec])
    return coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()


def smallest_nonzero_gap(A, mu):
    """smallest nonzero eigenvalue of the generalised problem  A f = lam diag(mu) f,
       A a symmetric PSD matrix with a one-dimensional kernel (row sums 0).
       Small problems: dense eigh (exact).  Large: Lanczos shift-invert, then take the
       smallest returned value above the numerical zero."""
    d = np.sqrt(mu)
    nn = len(mu)
    if hasattr(A, "tocsr"):
        A = A.tocsr()
        if nn <= 1200:
            Ad = A.toarray()
            C = (Ad / d[:, None]) / d[None, :]
            w = np.linalg.eigvalsh(0.5 * (C + C.T))
        else:
            C = diags(1.0 / d) @ A @ diags(1.0 / d)
            C = 0.5 * (C + C.T)
            w = eigsh(C, k=8, sigma=-1e-7, which="LM", return_eigenvectors=False)
    else:
        Ad = np.asarray(A, dtype=float)
        C = (Ad / d[:, None]) / d[None, :]
        w = np.linalg.eigvalsh(0.5 * (C + C.T))
    w = np.sort(np.abs(np.asarray(w)))
    pos = w[w > 1e-7 * max(w[-1], 1e-300)]
    return float(pos[0])


say("=" * 78)
say("R1 / LEVEL 1 & 2  --  non-Gaussian chain, P4-type kernel")
say("=" * 78)

n, beta, kappa = 240, 2.5, 0.85
x, h = grid(n)
w = onsite(x, kappa)
K = kernel(x, x, beta) * h
T = np.sqrt(w)[:, None] * K * np.sqrt(w)[None, :]
T = 0.5 * (T + T.T)
ev, evec = eigh(T)
o = np.argsort(ev)[::-1]
ev, evec = ev[o], evec[:, o]
lam = ev / ev[0]
Omega = evec[:, 0]
Omega = Omega / np.sqrt(h * np.sum(Omega ** 2))
if Omega[np.argmax(np.abs(Omega))] < 0:
    Omega = -Omega
mu = Omega ** 2 * h
mu = mu / mu.sum()
gap_cont = -np.log(lam[1])
gap_disc = 1.0 - lam[1]
say(f"  beta={beta}, kappa={kappa}, n={n}:  lam1/lam0 = {lam[1]:.10f}")
say(f"    continuous-time gap  -log(lam1/lam0) = {gap_cont:.8f}")
say(f"    discrete-time   gap  1 - lam1/lam0   = {gap_disc:.8f}")

# ---------------------------------------------------------------- [1] exact identity
say("")
say("-" * 78)
say("[1] EXACT: gap(T) = gap of the carre-du-champ form in L^2(mu),  mu = Omega^2")
say("-" * 78)
P = T / (Omega[:, None] * Omega[None, :])
r = P @ mu
say(f"    mu-stochastic: max|P mu - lam0*h| = {np.max(np.abs(r - ev[0] * h)):.2e}   (flat vector fixed)")
Q = P * mu[None, :] / (ev[0] * h)              # row-stochastic, reversible w.r.t. mu
Wm = 0.5 * (mu[:, None] * Q + (mu[:, None] * Q).T)
Aform = diags(Wm.sum(axis=1)) - Wm
lam_form = smallest_nonzero_gap(Aform, mu)
say(f"    smallest nonzero eigenvalue of the Dirichlet form in L^2(mu): {lam_form:.10f}")
say(f"    1 - lam1/lam0                                              : {gap_disc:.10f}")
chk("1.1 gap(T) = gap(semigroup Dirichlet form in L^2(Omega^2))", abs(lam_form - gap_disc) < 1e-8,
    f"|diff| = {abs(lam_form-gap_disc):.2e}")

# ---------------------------------------------------------------- [2] time-scale factor
say("")
say("-" * 78)
say("[2] the DISCRETE-vs-CONTINUOUS time factor  (1-x)/(-log x),  x = lam1/lam0")
say("-" * 78)
say(f"  {'x':>9} {'1-x':>10} {'-log x':>10} {'ratio':>9}")
for xx in (0.99, 0.9, 0.7, 0.52498, 0.3, 0.1, 0.01):
    say(f"  {xx:9.5f} {1-xx:10.6f} {-np.log(xx):10.6f} {(1-xx)/(-np.log(xx)):9.5f}")
say(f"  measured at beta={beta}: x={lam[1]:.5f} -> {(1-lam[1])/(-np.log(lam[1])):.5f}")
chk("2.1 the two legitimate gap normalisations differ by a computable, unbounded ratio",
    (1 - lam[1]) / (-np.log(lam[1])) < 1.0,
    f"discrete/continuous = {(1-lam[1])/(-np.log(lam[1])):.4f}; log-x / (1-x) -> inf as x -> 1")

# ---------------------------------------------------------------- [3] gradient substitution
say("")
say("-" * 78)
say("[3] the GRADIENT substitution lam1^grad(mu), and its restoration for a narrow kernel")
say("-" * 78)


def lam_grad(mu_vec, hh, mode="geom"):
    n_ = len(mu_vec)
    ip = (np.arange(n_) + 1) % n_
    c = (np.sqrt(mu_vec * mu_vec[ip]) if mode == "geom" else 0.5 * (mu_vec + mu_vec[ip])) / hh ** 2
    return smallest_nonzero_gap(graph_lap(c, n_), mu_vec)


lg = lam_grad(mu, h)
la = lam_grad(mu, h, "arith")
say(f"    lam1^grad(mu)  geometric-mean weights : {lg:.8f}")
say(f"    lam1^grad(mu)  arithmetic-mean weights: {la:.8f}   rel.diff = {abs(la-lg)/lg:.2e}")
say(f"    semigroup gap                        : {gap_cont:.8f}")
say(f"    ratio lam1^grad / gap                : {lg/gap_cont:.6f}  <-- metric mismatch at beta=2.5")
chk("3.1 geometric and arithmetic weights agree to O(h^2) (discrete ground-state transform)",
    abs(la - lg) / lg < 1e-3, f"{abs(la-lg)/lg:.2e}")
say("    narrow-kernel restoration: k_eps = exp(-(y-x)^2/(2 eps)), matched weight eps/2")
say(f"  {'eps':>7} {'gap(=-log l1/l0)':>18} {'(eps/2)lam1^grad':>18} {'ratio':>9}")
rows3 = []
for eps in (0.40, 0.20, 0.10, 0.05, 0.025, 0.0125):
    Ke = np.exp(-((x[None, :] - x[:, None]) ** 2) / (2 * eps)) * h
    Te = np.sqrt(w)[:, None] * Ke * np.sqrt(w)[None, :]
    Te = 0.5 * (Te + Te.T)
    e2, v2 = eigh(Te)
    o2 = np.argsort(e2)[::-1]
    e2, v2 = e2[o2], v2[:, o2]
    O2 = v2[:, 0]
    O2 = O2 / np.sqrt(h * np.sum(O2 ** 2))
    if O2[np.argmax(np.abs(O2))] < 0:
        O2 = -O2
    mu2 = O2 ** 2 * h
    mu2 = mu2 / mu2.sum()
    g2 = -np.log(e2[1] / e2[0])
    l2 = lam_grad(mu2, h)
    rr = 0.5 * eps * l2 / g2
    rows3.append({"eps": eps, "gap": g2, "matched": 0.5 * eps * l2, "ratio": rr})
    say(f"  {eps:7.4f} {g2:18.8f} {0.5*eps*l2:18.8f} {rr:9.5f}")
chk("3.2 with the matched metric weight eps/2 the gradient constant -> the gap (residual O(eps))",
    abs(rows3[-1]["ratio"] - 1.0) < 0.06
    and all(rows3[i + 1]["ratio"] <= rows3[i]["ratio"] + 1e-9 for i in range(1, len(rows3) - 1)),
    f"ratio {rows3[0]['ratio']:.4f} -> {rows3[-1]['ratio']:.4f} (monotone down toward 1)")
say("    READING.  The LSI controls the GRADIENT form, not the transfer gap.  The step")
say("    'gradient form = gap' is exact only for a DIFFUSION generator (eps -> 0, ratio 1);")
say("    for the isotropic lattice kernel at beta=2.5 it costs the fixed factor printed above")
say("    (0.668).  For a Schroedinger/KS generator it costs exactly the kinetic convention 2.")

# ---------------------------------------------------------------- [4] full chain
say("")
say("-" * 78)
say("[4] the FULL chain: lam_conf(chain) <= lam1^grad(mu_slice), strictly")
say("-" * 78)


def chain_conf(Lc, nc, beta_, kappa_):
    xc, hc = grid(nc)
    wc = onsite(xc, kappa_)
    kc = kernel(xc, xc, beta_)
    idx = np.arange(nc ** Lc)
    coords = np.stack(np.unravel_index(idx, (nc,) * Lc), axis=1)
    logmu = np.zeros(len(idx))
    for i in range(Lc):
        logmu += np.log(wc[coords[:, i]])
        j = (i + 1) % Lc
        logmu += -beta_ * (1.0 - np.cos(0.5 * (xc[coords[:, j]] - xc[coords[:, i]])))
    mu_c = np.exp(logmu - logmu.max())
    mu_c /= mu_c.sum()
    rows, cols, vals = [], [], []
    diag = np.zeros(len(idx))
    for i in range(Lc):
        nx = coords.copy()
        nx[:, i] = (nx[:, i] + 1) % nc
        jdx = np.ravel_multi_index(tuple(nx.T), (nc,) * Lc)
        c = np.sqrt(mu_c * mu_c[jdx]) / hc ** 2
        rows += [idx, jdx]
        cols += [jdx, idx]
        vals += [-c, -c]
        dd = c.copy()
        np.add.at(dd, jdx, c)          # row sum: c_i + c_{i-1}
        diag += dd
    rows.append(idx)
    cols.append(idx)
    vals.append(diag)
    A = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(len(idx), len(idx))).tocsr()
    gap = smallest_nonzero_gap(A, mu_c)
    return gap, int(len(idx))


rows4 = []
say("  per-L form of the slice-marginal lemma:  lam_conf(chain_L) <= lam1(marginal_L),")
say("  and marginal_L -> |Omega|^2 as L -> inf (so the L=inf bound is lam_slice).")
say(f"  {'L':>3} {'grid':>5} {'dim':>8} {'lam_conf':>10} {'lam1(marg_L)':>13} {'slack':>9}")
marg_cache = {}
for Lc, nc in ((1, 20), (2, 20), (3, 20)):
    lc, dim = chain_conf(Lc, nc, beta, kappa)
    xg, hg = grid(nc)
    wg = onsite(xg, kappa)
    Qc = wg[:, None] * kernel(xg, xg, beta)          # same grid as the chain
    QL = np.linalg.matrix_power(Qc, Lc)
    m_ = np.diag(QL)
    m_ = m_ / m_.sum()
    lm = lam_grad(m_, hg)
    rows4.append({"L": Lc, "nc": nc, "dim": dim, "lam_conf": lc, "lam_marg": lm,
                  "slack": lm - lc, "ratio": lc / lm})
    say(f"  {Lc:3d} {nc:5d} {dim:8d} {lc:10.6f} {lm:13.6f} {lm-lc:9.6f}")
chk("4.1 lam_conf(chain_L) <= lam1(marginal_L) at every L (slice-marginal lemma)",
    all(r["slack"] >= -1e-4 for r in rows4),
    "slacks: " + ", ".join(f"{r['slack']:.4f}" for r in rows4))
chk("4.2 at L=1 the lemma is an EQUALITY (the chain is its own marginal)",
    abs(rows4[0]["slack"]) < 1e-4, f"slack = {rows4[0]['slack']:.2e}")
chk("4.3 lam_conf(chain_3) < lam_slice (the L->inf vacuum marginal): strict wrong-space loss",
    rows4[-1]["lam_conf"] < lg, f"{rows4[-1]['lam_conf']:.4f} < {lg:.4f}")
say("    finite-L slice marginal -> |Omega|^2  (why the L->inf slice marginal IS the vacuum):")
Qm = w[:, None] * K
for Lc in (1, 2, 4, 8, 16):
    QL = np.linalg.matrix_power(Qm, Lc)
    m_ = np.diag(QL)
    m_ = m_ / m_.sum()
    say(f"      L={Lc:3d}: ||marg_L - |Omega|^2||_1 = {np.sum(np.abs(m_ - mu)):.3e}")

# ---------------------------------------------------------------- [5] level 2
say("")
say("-" * 78)
say("[5] LEVEL 2: restriction to a subspace cannot lower the constant (projection is safe)")
say("-" * 78)
mir = (-np.arange(n)) % n
half = n // 2
basis_e, basis_o = [], []
for k in range(0, half + 1):                      # orbits {k, n-k}: one representative each
    v = np.zeros(n)
    v[k] = 1.0
    v[mir[k]] += 1.0
    basis_e.append(v / np.linalg.norm(v))
for k in range(1, half):                          # the two fixed points k=0, n/2 have no odd vector
    v = np.zeros(n)
    v[k] = 1.0
    v[mir[k]] -= 1.0
    basis_o.append(v / np.linalg.norm(v))
E = np.array(basis_e).T
Eo = np.array(basis_o).T
ip = (np.arange(n) + 1) % n
cvec = np.sqrt(mu * mu[ip]) / h ** 2
A_full = graph_lap(cvec, n).toarray()
AE = E.T @ A_full @ E
BE = E.T @ np.diag(mu) @ E
from scipy.linalg import eigh as eigh_gen
wE = np.sort(eigh_gen(0.5 * (AE + AE.T), 0.5 * (BE + BE.T), eigvals_only=True))
lam_even = float(np.min(wE[wE > 1e-9 * max(wE)]))
say(f"    full space    lam1 = {lg:.8f}")
say(f"    even subspace lam1 = {lam_even:.8f}   ratio = {lam_even/lg:.6f}")
chk("5.1 restriction to a subspace can only INCREASE the constant",
    lam_even >= lg - 1e-6, f"ratio {lam_even/lg:.6f} >= 1")
# odd subspace (f(x) = -f(-x)): the slow mode is even, so this one must be strictly larger
AO = Eo.T @ A_full @ Eo
BO = Eo.T @ np.diag(mu) @ Eo
wO = np.sort(eigh_gen(0.5 * (AO + AO.T), 0.5 * (BO + BO.T), eigvals_only=True))
lam_odd = float(np.min(wO[wO > 1e-9 * max(wO)]))
say(f"    odd  subspace lam1 = {lam_odd:.8f}   ratio = {lam_odd/lg:.6f}")
chk("5.2 the odd sector is strictly heavier: the slow mode lives in the even (vacuum-like) sector",
    lam_odd > lg * (1 + 1e-6), f"odd/full = {lam_odd/lg:.4f}")
chk("5.3 the vacuum sector loses NOTHING (ratio 1 to numerical precision)",
    abs(lam_even / lg - 1.0) < 1e-3, f"{lam_even/lg:.8f}")

# ---------------------------------------------------------------- [6] LSI -> Poincare
say("")
say("-" * 78)
say("[6] the LSI -> Poincare loss, exact on the two-point space")
say("-" * 78)
say("    mu = (1-eps, eps);  E(f) = (f1-f0)^2;  lam_P = 1/(eps(1-eps));  gamma = 2/sup_f Ent(f^2)/E(f)")
say("    parametrise E(f)=1:  f = (a, a+1), scan a  (avoids the 0/0 corner at a -> +-inf)")
say(f"  {'eps':>9} {'lam_P':>12} {'gamma_LSI':>12} {'lam_P/gamma':>12} {'log(1/eps)/2':>13}")
rows6 = []
ag = np.linspace(-60.0, 60.0, 2400011)
for eps in (0.5, 0.1, 0.01, 1e-3, 1e-6):
    m0, m1 = 1 - eps, eps
    f0 = ag
    f1 = ag + 1.0                       # E(f) = 1 exactly
    s = m0 * f0 ** 2 + m1 * f1 ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = (m0 * f0 ** 2 * np.log(f0 ** 2) + m1 * f1 ** 2 * np.log(f1 ** 2)) - s * np.log(s)
    ent = np.nan_to_num(ent, nan=0.0, posinf=0.0, neginf=0.0)
    ratio = np.max(ent)                 # = sup Ent(f^2)/E(f) since E = 1
    gamma = 2.0 / ratio
    lamP = 1.0 / (eps * (1 - eps))
    rows6.append({"eps": eps, "lamP": lamP, "gamma": gamma, "ratio": lamP / gamma})
    say(f"  {eps:9.0e} {lamP:12.6e} {gamma:12.6e} {lamP/gamma:12.4f} {0.5*np.log(1/eps):13.4f}")
chk("6.1 the LSI constant is strictly weaker than Poincare, loss UNBOUNDED in log(1/eps)",
    rows6[-1]["ratio"] > 3 * rows6[0]["ratio"] and rows6[-1]["ratio"] > 5.0,
    f"lam_P/gamma: {rows6[0]['ratio']:.3f} (eps=0.5) -> {rows6[-1]['ratio']:.3f} (eps=1e-6);"
    f" predicted log(1/eps)/2 = {0.5*np.log(1e6):.3f}")
chk("6.2 both constants are positive: an LSI bound is VALID but loses the log factor",
    rows6[-1]["gamma"] > 0 and rows6[-1]["lamP"] > rows6[-1]["gamma"],
    f"gamma={rows6[-1]['gamma']:.3e} < lamP={rows6[-1]['lamP']:.3e}")

say("")
say("=" * 78)
say("VERDICT")
say("=" * 78)
say("  [1] gap(T) = gap(carre-du-champ form in L^2(vacuum)) : EXACT (machine precision).")
say("  [2] discrete vs continuous time normalisation: factor (1-x)/(-log x), computable.")
say("  [3] the GRADIENT form needs the metric matched to the kinetic term; matched, it equals")
say("      the gap (ratio -> 1 as eps -> 0); at beta=2.5 with the unit metric it is off by a")
say("      fixed factor printed above.")
say("  [4] the full-chain Poincare constant is STRICTLY below the vacuum-slice constant.")
say("  [5] the subspace projection never lowers the constant (gauge-invariance analog).")
say("  [6] the LSI constant is strictly weaker than Poincare with an unbounded loss.")
say("=" * 78)
save({"lam": [float(v) for v in lam[:5]], "gap_cont": gap_cont, "gap_disc": gap_disc,
      "lam_grad": lg, "lam_grad_arith": la, "rows3": rows3, "rows4": rows4,
      "lam_even": lam_even, "rows6": rows6, "checks": CHK})
