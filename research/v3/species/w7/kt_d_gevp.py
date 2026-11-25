"""kt_d_gevp.py -- KT-RHO2 T1 (second method + tightness) and T2 (multiplicity), plus the
read-only 4D ledger.

D1  two independent methods for rho^2: direct diagonalisation vs GEVP (Ritz values of a
    coarse trial space), with the error as a function of the trial dimension.
D2  tightness of  rho^2 <= t/e1 - 1 ,  t = T2(V) = tr(G^-1 M)   (RB's rank-2 trace).
    How loose is the constructible bound, and what floor does `e1 >= t/(1+C^2)` demand?
D3  bottom multiplicity / near-degeneracy pattern table, assembled from kt_a/b/c JSONs.
D4  read-only 4D ledger: what the stored w7 data actually contain.
Convention: e1 = lowest generator eigenvalue of the tower (vacuum already removed when the
model has one); rho^2 = e2/e1.
"""
import json, glob, os, numpy as np
from scipy.linalg import eigh_tridiagonal, eigh

JSON, LOG = "kt_d_gevp.json", []
OUT = {}


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.append(s)


def flush():
    with open(JSON, "w", encoding="utf-8") as f:
        json.dump(OUT, f, indent=1)
    with open("kt_d_gevp.log", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


def chain(L):
    return np.full(L, 2.0), np.full(L - 1, -1.0)


def gevp_ritz(A, basis):
    """Ritz values of  M x = nu G x  on span(basis),  basis: n x m."""
    B = np.asarray(basis, float)
    G = B.T @ B
    M = B.T @ (A @ B)
    return np.sort(eigh(M, G, eigvals_only=True))


P("=" * 100)
P("D1  method A (direct) vs method B (GEVP): rho^2 on a path chain")
P("=" * 100)
d1 = {}
for L in [40, 200]:
    d, e = chain(L)
    A = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
    lam = np.sort(eigh_tridiagonal(d, e, eigvals_only=True))
    mu = lam                                    # Dirichlet-ended: no vacuum in the tower
    rho_direct = mu[1] / mu[0]
    P(f"  L={L}: direct e1={mu[0]:.10e} e2={mu[1]:.10e} e3={mu[2]:.10e}"
      f"  rho^2={rho_direct:.10f}")
    rows = {}
    for m in [2, 4, 8, 16, 32]:
        blocks = np.array_split(np.arange(L), m)
        basis = np.zeros((L, m))
        for k, bl in enumerate(blocks):
            basis[bl, k] = 1.0
        nu = gevp_ritz(A, basis)
        r_gevp = nu[1] / nu[0]
        rows[m] = dict(rho2=float(r_gevp), err=float(abs(r_gevp - rho_direct)))
        P(f"    coarse blocks m={m:<3d}: nu1={nu[0]:.6e} nu2={nu[1]:.6e}"
          f"  rho^2_GEVP={r_gevp:.8f}  |err|={abs(r_gevp-rho_direct):.2e}")
    W = np.linalg.eigh(A)[1][:, 0:3]
    nu = gevp_ritz(A, W)
    P(f"    exact low space : rho^2_GEVP={nu[1]/nu[0]:.10f}"
      f"  err={abs(nu[1]/nu[0]-rho_direct):.2e}")
    d1[f"L{L}"] = dict(direct=float(rho_direct), gevp=rows, exact_space=float(nu[1] / nu[0]))
OUT["D1_two_methods"] = d1
flush()

P("=" * 100)
P("D2  tightness:  t = T2(V) = tr(G^-1 M)  vs  the optimal t* = e1+e2")
P("    bound rho^2 <= t/e1 - 1 ;  RB-2'' threshold  e1 >= t/(1+C^2),  C=2")
P("=" * 100)


def tightness(A, basis, mu, tag):
    nu = gevp_ritz(A, basis)
    t = float(np.sum(nu))
    t_star = mu[0] + mu[1]
    rho = mu[1] / mu[0]
    bound = t / mu[0] - 1.0
    thr = t / 5.0
    d = dict(t=t, t_star=float(t_star), slack=float(t / t_star), rho2=float(rho),
             bound=float(bound), loose=float(bound / rho), thr=float(thr),
             thr_over_e1=float(thr / mu[0]))
    P(f"  {tag:<30s} t={t:.6e} t*={t_star:.6e} t/t*={t/t_star:9.3f}"
      f"  bound={bound:10.3f} (true {rho:.4f}, x{bound/rho:6.1f} loose)"
      f"  t/5={thr/mu[0]:7.2f} x e1")
    return d


d2 = {}
L = 40
d, e = chain(L)
A = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
mu = np.sort(eigh_tridiagonal(d, e, eigvals_only=True))
for b in [1, 2, 4, 8, 20]:
    basis = np.zeros((L, 2))
    basis[max(0, L // 2 - b):L // 2, 0] = 1
    basis[L // 2:L // 2 + b, 1] = 1
    basis /= np.linalg.norm(basis, axis=0)
    d2[f"adjacent_blocks_b{b}"] = tightness(A, basis, mu, f"adjacent blocks b={b}")
for r in [1.0, 2.0, 4.0]:
    j = np.arange(L)
    o1 = np.exp(-0.5 * ((j - L / 2) / r) ** 2)
    o2 = np.exp(-0.5 * ((j - L / 2) / (2 * r)) ** 2)
    basis = np.stack([o1, o2], 1)
    basis /= np.linalg.norm(basis, axis=0)
    d2[f"two_smear_r{r:g}"] = tightness(A, basis, mu, f"two smearing radii r={r:g}")
W = np.linalg.eigh(A)[1][:, 0:2]
d2["exact_low_space"] = tightness(A, W, mu, "exact low space (optimal)")
OUT["D2_tightness"] = d2
flush()

P("  double well Lam=20 (N=4000): the two-well coarse space and the tunnelling mode")
N = 4000
x = np.linspace(-3, 3, N)
h = x[1] - x[0]
V = 20.0 * (x ** 2 - 1) ** 2
dW = 2.0 / h ** 2 + V
eW = np.full(N - 1, -1.0 / h ** 2)
ev = eigh_tridiagonal(dW, eW, eigvals_only=True, select="i", select_range=(0, 5))
muw = ev - ev[0]
P(f"    mu1={muw[1]:.4e} mu2={muw[2]:.6f} rho^2={muw[2]/muw[1]:.3f}")
AW = np.diag(dW) + np.diag(eW, 1) + np.diag(eW, -1)
basis = np.zeros((N, 2))
basis[np.abs(x + 1) < 0.4, 0] = 1
basis[np.abs(x - 1) < 0.4, 1] = 1
basis /= np.linalg.norm(basis, axis=0)
d2["double_well_two_wells"] = tightness(AW, basis, muw[1:], "double well, two-well basis")
OUT["D2_tightness"] = d2
flush()

P("=" * 100)
P("D3  bottom multiplicity / near-degeneracy pattern (T2)")
P("=" * 100)
tab = {}
for f, key in [("kt_a_exact.json", "A1_chains"), ("kt_a_exact.json", "A2_box_BC"),
               ("kt_a_exact.json", "A3_plaquette_2dYM"), ("kt_a_exact.json", "A5_rotor_U1"),
               ("kt_a_exact.json", "A6_radial"), ("kt_a_exact.json", "A7_discrete_star")]:
    if not os.path.exists(f):
        continue
    J = json.load(open(f, encoding="utf-8"))
    for k, v in J.get(key, {}).items():
        if v.get("rho2") is None:
            continue
        tab[f"{key}:{k}"] = dict(rho2=v["rho2"], mult=v.get("mult"), e1=v.get("e1"),
                                 e2=v.get("e2"), near=v.get("near_deg"))
if os.path.exists("kt_c_geom.json"):
    J = json.load(open("kt_c_geom.json", encoding="utf-8"))
    for k, v in J.get("C1_graphs", {}).items():
        tab[f"C1_graphs:{k}"] = dict(rho2=v["merged"], mult=v.get("mult"),
                                     split=v.get("split"), lam=v.get("lam"))
    for k, v in J.get("C2_domains", {}).items():
        tab[f"C2_domains:{k}"] = dict(rho2=v["rho2"], mult=None, split=None, lam=v.get("lam"))
    for k, v in J.get("C3_double_well", {}).items():
        tab[f"C3_double_well:{k}"] = dict(rho2=v["rho2"], mult=None,
                                          split=v.get("pair_split"))
for k, v in tab.items():
    sp = v.get("split")
    near = (v.get("near") is not None and v["near"] < 1e-2) or (sp is not None and sp < 1e-2)
    P(f"  {k:<44s} rho^2={v['rho2']:11.5f}  mult={str(v.get('mult')):<18s}"
      f"  near_degen={'YES' if near else 'no'}")
OUT["D3_pattern"] = tab
flush()

P("=" * 100)
P("D4  read-only 4D ledger: what the stored w7 data actually contain")
P("=" * 100)
led = {}
files = [f for f in glob.glob("*.json")] + [f for f in glob.glob("../w7/*.json")]
for f in sorted(set(files)):
    try:
        J = json.load(open(f, encoding="utf-8"))
    except Exception:
        continue
    ks = [str(k)[:56] for k in (list(J.keys())[:8] if isinstance(J, dict) else [])]
    flat = json.dumps(J)[:6000].lower()
    gevp_like = any(t in flat for t in ("gevp", "corr_mat", "cmat", "nmode", "eigvec"))
    led[os.path.basename(f)] = dict(keys=ks, multi_operator=bool(gevp_like))
for k, v in led.items():
    P(f"  {k:<32s} multi-op? {'YES' if v['multi_operator'] else 'no ':<4s} keys={v['keys']}")
OUT["D4_4d_ledger"] = led
flush()
P("done. json ->", JSON)
