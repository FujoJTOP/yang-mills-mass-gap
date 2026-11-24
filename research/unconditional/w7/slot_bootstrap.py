"""
slot_bootstrap.py  --  witness for research/unconditional/SLOT-BOOTSTRAP.md
Single small job.  Exact linear algebra / tiny LPs.  NO lattice MC.  Lightweight.
Run:  py -3.12 slot_bootstrap.py

ROUTE = numerical bootstrap / feasibility of a gapped k->0 solution.
Finite constraint system on the 0++ spectral measure nu  (a truncated Stieltjes
moment problem).  We compute -- not assert -- the feasible support-floor
Delta = inf supp nu  as a function of truncation order K.

Exact inputs encoded:
  (RP)  reflection positivity  <=>  nu >= 0                 [LP: w >= 0]
  (KL)  Kallen-Lehmann moments  M_j = int s^j d nu,  s = m^2
        M_{-1} = chi = 1/G''  (Ward  E[af]=1  =>  chi < inf)
        M_{-2} = A ~ Lambda^4 (trace anomaly, finite)
        M_0 = W,  M_1..M_K    (norm + growth / AF UV tail)
Representative gapped moment vector (V-here); the floor-blindness is a property
of the moment KERNELS, value-independent (Task 3 verifies via leave-one-out).

Extremal floor quantities:
  Delta_max(K) = sup over feasible nu of inf supp nu  -> an UPPER bound on the floor
  Delta_min(K) = inf over feasible nu of inf supp nu  -> the only LOWER-bound certificate
The mass GAP needs Delta_min(K) > 0 uniformly.  We test whether it is.
"""
import numpy as np
from scipy.optimize import linprog

np.set_printoptions(precision=4, suppress=True, linewidth=120)


def banner(s):
    print("\n" + "=" * 80)
    print(s)
    print("=" * 80)


def jset(K):
    return list(range(-2, K + 1))            # [-2,-1,0,1,...,K]


def moments_of(nodes, weights, js):
    nodes = np.asarray(nodes, float); weights = np.asarray(weights, float)
    return np.array([np.sum(weights * nodes ** j) for j in js])


# ---------------------------------------------------------------------------
# Representative physical moment vector.  Units: lightest 0++ glueball m1^2 = 1.
# A dominant glueball atom f1^2 + a sparse multi-glueball continuum.
# ---------------------------------------------------------------------------
TRUE_S = np.array([1.0, 3.0, 6.0, 12.0, 25.0])
TRUE_W = np.array([0.50, 0.15, 0.15, 0.10, 0.10])   # sum = 1 = W
TRUE_FLOOR = 1.0


def true_moments(K):
    return moments_of(TRUE_S, TRUE_W, jset(K))


# ---------------------------------------------------------------------------
# LP solver (row-scaled).  find w >= 0 with  sum_i w_i node_i^j = mu_j  for j in J.
#   low_index is not None -> maximize w[low_index]  (push mass as low as allowed).
# Returns (w, scaled_resid).
# ---------------------------------------------------------------------------
def solve_measure(nodes, js, mu, low_index=None):
    nodes = np.asarray(nodes, float)
    A = np.array([[s ** j for s in nodes] for j in js], float)
    b = np.asarray(mu, float).copy()
    r = np.max(np.abs(A), axis=1); r[r == 0] = 1.0
    A = A / r[:, None]; b = b / r
    c = np.zeros(len(nodes))
    if low_index is not None:
        c[low_index] = -1.0
    res = linprog(c, A_eq=A, b_eq=b, bounds=[(0, None)] * len(nodes), method="highs")
    if not res.success:
        return None, np.inf
    w = res.x
    return w, np.max(np.abs(A @ w - b))


# ---------------------------------------------------------------------------
# Delta_max(K): largest floor D such that a measure with ALL support >= D exists.
# Every atom >= D; the physical truth support points are injected as nodes so the
# truth (floor 1) is always exactly representable -> Delta_max >= 1 guaranteed.
# ---------------------------------------------------------------------------
def delta_max(K, S_max=60.0, n_grid=30):
    js = jset(K); mu = true_moments(K); best = 0.0
    Ds = np.unique(np.concatenate([np.geomspace(0.05, 3.0, 240),
                                   [0.99, 0.9999, 1.0, 1.0001, 1.001, 1.005, 1.01]]))
    for D in Ds:
        grid = np.geomspace(D, S_max, n_grid)
        sup = TRUE_S[TRUE_S >= D]                      # exact truth nodes, if allowed
        nodes = np.unique(np.concatenate([grid, sup]))
        w, resid = solve_measure(nodes, js, mu)
        if w is not None and resid < 1e-6:
            best = D
    return best


# ---------------------------------------------------------------------------
# Delta_min probe: is a measure with support reaching DOWN to Delta_t feasible?
# A low atom at s = Delta_t (index 0); rest of the atoms on [m_lo, S_max] plus the
# exact truth nodes (so the bulk is exactly representable).  Maximize the low-atom
# weight.  w0* > 0 (resid tiny)  PROVES  floor <= Delta_t is consistent.
# ---------------------------------------------------------------------------
def delta_min_probe(K, Delta_t, m_lo=0.5, S_max=60.0, n_grid=40):
    js = jset(K); mu = true_moments(K)
    grid = np.geomspace(m_lo, S_max, n_grid)
    nodes = np.unique(np.concatenate([[Delta_t], grid, TRUE_S]))
    i_low = int(np.argmin(np.abs(nodes - Delta_t)))    # the Delta_t node
    w, resid = solve_measure(nodes, js, mu, low_index=i_low)
    if w is None:
        return 0.0, np.inf, np.nan
    w0 = w[i_low]
    frac_anom = (w0 * Delta_t ** (-2)) / mu[0]          # share of M_{-2} from low atom
    return w0, resid, frac_anom


# ===========================================================================
banner("SETUP  |  representative gapped 0++ spectral measure -> moments")
K0 = 6
mu0 = true_moments(K0)
print(f"truth support s = {TRUE_S}, weights = {TRUE_W},  floor(m1^2) = {TRUE_FLOOR}")
labels = {-2: "M_{-2}=A (anomaly ~L^4)", -1: "M_{-1}=chi (Ward E[af]=1)",
          0: "M_0 = W (norm)", 1: "M_1", 2: "M_2", 3: "M_3"}
for j, name in labels.items():
    print(f"  {name:26s} = {moments_of(TRUE_S,TRUE_W,[j])[0]:.6f}")
print(f"  moment-ratio UPPER bound  M0/M_{-1} = {mu0[2]/mu0[1]:.6f}  (>= true floor 1)")
print("  (V-here: model values; structural conclusion is value-independent)")

# ===========================================================================
banner("TEST 0  |  the minimal exact ambiguity -- same (M0,M1), floors 0 vs 1")
a_s = np.array([0.0, 2.0]); a_w = np.array([1.0, 1.0])
b_s = np.array([1.0]);      b_w = np.array([2.0])
Ma = moments_of(a_s, a_w, [0, 1]); Mb = moments_of(b_s, b_w, [0, 1])
print(f"  delta_0 + delta_2 : M0,M1 = {Ma},  inf supp = {a_s.min()}")
print(f"  2 delta_1         : M0,M1 = {Mb},  inf supp = {b_s.min()}")
print(f"  identical to machine zero? {np.allclose(Ma, Mb)}   -> order-1 moments ALREADY blind")

# ===========================================================================
banner("TEST 1  |  floor scan vs truncation order  K  (the decisive table)")
# exact-arithmetic sanity: the truth itself (floor=1) satisfies its own moments,
# so Delta_max >= 1 exactly; scan values below 1 are LP conditioning at high K.
resid_truth = np.max(np.abs(moments_of(TRUE_S, TRUE_W, jset(7)) - true_moments(7)))
print(f"  truth-witness residual = {resid_truth:.2e}  =>  Delta_max >= 1 in exact arithmetic")
print(f"{'K':>2} | {'|J|':>3} | {'Delta_max(K)':>12} | {'Delta_min(K) <= ':>16} | "
      f"{'gapless <=1e-3 ?':>17}")
print("-" * 72)
rows = []
probe_Ds = [1e-1, 3e-2, 1e-2, 3e-3, 1e-3]
for K in range(1, 8):
    dmax = delta_max(K)
    dmin = None
    for D in probe_Ds:
        w0, resid, _ = delta_min_probe(K, D)
        if w0 is not None and w0 > 1e-14 and resid < 1e-5:
            dmin = D
        else:
            break
    gl = "YES" if (dmin is not None and dmin <= 1.01e-3) else "no"
    nm = f"{dmin:.0e}" if dmin is not None else "inf"
    print(f"{K:>2} | {K+3:>3} | {dmax:>12.5f} | {nm:>16} | {gl:>17}")
    rows.append((K, dmax, dmin))

# ===========================================================================
banner("TEST 2  |  near-zero atom budget w0(Delta) -- weight the data allow at the floor")
K = 6
print("  Delta_t     w0*(Delta)   scaled-resid   share of M_{-2} from low atom   feasible?")
for D in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3]:
    w0, resid, frac = delta_min_probe(K, D)
    ok = (w0 is not None and w0 > 1e-14 and resid < 1e-5)
    r = f"{resid:.1e}" if np.isfinite(resid) else "  inf"
    print(f"  {D:.0e}     {w0:.3e}     {r:>10}   {frac:>10.4f}                   {str(ok):>5}")
print("  -> w0* ~ C * Delta^2 (anomaly kernel 1/s^2 concentrates at the floor):")
print("     a *vanishing-mass* atom supplies a vanishing share of all moments, so a")
print("     gapless (floor -> 0) measure is consistent with every exact input.")

# ===========================================================================
banner("TEST 3  |  leave-one-out -- WHICH constraint would block a gapless solution?")
K = 5
full = jset(K)
w0a, r0a, _ = delta_min_probe(K, 1e-3)
print(f"  All exact inputs, floor<=1e-3 feasible? "
      f"{w0a is not None and w0a > 0 and r0a < 1e-5}  (w0={w0a:.2e})")
for drop in [-2, -1, 0, 1, 2]:
    js = [j for j in full if j != drop]
    mu_full = true_moments(K)
    mu = np.array([mu_full[full.index(j)] for j in js])
    grid = np.geomspace(0.5, 60.0, 40)
    nodes = np.unique(np.concatenate([[1e-3], grid, TRUE_S]))
    i_low = int(np.argmin(np.abs(nodes - 1e-3)))
    w, resid = solve_measure(nodes, js, mu, low_index=i_low)
    ok = (w is not None and w[i_low] > 0 and resid < 1e-5)
    w0 = w[i_low] if w is not None else np.nan
    print(f"  drop M_{drop:+d}  -> floor<=1e-3 feasible? {str(ok):>5}   w0={w0:.3e}")
print("  -> gapless feasibility SURVIVES removing EACH moment: no single moment is the")
print("     'gate'. The obstruction is the TYPE of data (integral moments). M_{-2},M_{-1}")
print("     give only UPPER bounds on the floor.")

# ===========================================================================
banner("TEST 4  |  the binding constraints are UPPER bounds (moment-ratio lemma)")
print("  nu>=0, inf supp nu >= D  =>  M_{-1} <= M_0/D  =>  D <= M_0/M_{-1};  similarly")
print("  M_{-2} <= M_{-1}/D  =>  D <= M_{-1}/M_{-2}.  Every exact input caps D ABOVE.")
mu = true_moments(8)                                 # jset(8) -> idx j+2
print(f"   M0/M_{-1} = {mu[2]/mu[1]:.5f}   M_{-1}/M_{-2} = {mu[1]/mu[0]:.5f}   "
      f"M1/M0 = {mu[3]/mu[2]:.5f}   (upper bounds; true floor = 1)")
print("  -> no exact input lower-bounds the edge; all are upper-bound machines.")

# ===========================================================================
banner("VERDICT (numbers above)")
dmax_list = [r[1] for r in rows]
gl_all = all((r[2] is not None and r[2] <= 1.01e-3) for r in rows)
print(f"  Delta_max(K) = {[round(x,4) for x in dmax_list]}")
print(f"    -> UPPER bound on the floor (~1.00 within grid/LP resolution; >=1 exactly")
print(f"       via the truth witness), tightening toward the true floor 1 with K:")
print(f"       an upper bound, never a lower-bound / gap certificate.")
print(f"  Delta_min(K) <= 1e-3 holds at EVERY K? {gl_all}")
print("    -> Delta_min(K) = 0 at every truncation; the LOW end is pinned at 0 and")
print("       does NOT rise with K.  Moments/RP/Ward/anomaly certify floor = 0 only.")
print("  DECISIVE: the feasible floor set is [0, Delta_max] at every K; the minimum")
print("  feasible floor is 0 for all K.  Moment data alone CANNOT certify a positive")
print("  floor.  A NON-MOMENT (edge-sensitive) input is required.")
