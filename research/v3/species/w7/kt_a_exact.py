"""kt_a_exact.py -- KT-RHO2 T1/T2: rho^2 = e2/e1 on exactly solvable towers.

Little Code Sauce / Round 3 witness. Single process, deterministic, numpy+scipy.
Writes kt_a_exact.json after every model group (crash-safe) and kt_a_exact.log.

Purpose: fire-test C=2 (`rho^2 <= 4`) on EXACT towers, including the boundary-condition
families (Dirichlet/Neumann/periodic/mixed) which decide the constant.
Labels: [T] analytic spectrum, [V] numerical diagonalisation.
"""
import json, numpy as np
from scipy.linalg import eigh_tridiagonal, eigh
from scipy.special import jn_zeros

OUT, LOG = {}, []
JSON = "kt_a_exact.json"
LOG = []


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.append(s)


def flush():
    with open(JSON, "w", encoding="utf-8") as f:
        json.dump(OUT, f, indent=1)
    with open("kt_a_exact.log", "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


def group(mu, tol=1e-8):
    """merge (numerically) equal eigenvalues -> (distinct values, multiplicities)"""
    vals, mult = [], []
    for v in np.sort(np.asarray(mu, float)):
        if vals and abs(v - vals[-1]) <= tol * max(1.0, abs(v)):
            mult[-1] += 1
        else:
            vals.append(v)
            mult.append(1)
    return np.array(vals), np.array(mult)


def stats(mu, tag, vacuum=True, tol=1e-8, note="", merge=True):
    """rho^2 = e2/e1 with e_n = distinct nonzero (or lowest) generator eigenvalues."""
    raw = np.sort(np.asarray(mu, float))
    vals, mult = group(raw, tol=tol)
    if vacuum:                      # a zero mode is the vacuum, not part of the tower
        keep = vals > 1e-10 * max(1.0, abs(vals).max() if len(vals) else 1.0)
        vals, mult = vals[keep], mult[keep]
    if not merge:                   # raw = keep every (possibly degenerate) level
        vals, mult = raw[raw > 1e-10 * max(1.0, abs(raw).max())], np.ones(len(raw))
    if len(vals) < 2:
        P(f"  {tag:<34s} only {len(vals)} distinct level(s) -> rho^2 UNDEFINED  {note}")
        return dict(tag=tag, note=note, distinct=[float(v) for v in vals],
                    mult=[int(m) for m in mult], e1=float(vals[0]) if len(vals) else None,
                    e2=None, e3=None, rho2=None, rho2_raw=None, r3=None,
                    near_deg=None, Cmin=None, window_C2=None)
    e1, e2, e3 = vals[0], vals[1], (vals[2] if len(vals) > 2 else np.nan)
    rho2 = e2 / e1
    d = dict(
        tag=tag, note=note,
        raw_first=list(np.round(raw[:6], 10)),
        distinct=[float(v) for v in vals[:5]],
        mult=[int(m) for m in mult[:5]],
        e1=float(e1), e2=float(e2), e3=float(e3) if e3 == e3 else None,
        rho2=float(rho2),
        rho2_raw=float(raw[1] / raw[0]) if (vacuum and len(raw) > 1 and raw[0] > 0) else None,
        r3=float(e3 / e2) if e3 == e3 else None,
        near_deg=float(abs(e2 - e1) / e1),
        Cmin=float(np.sqrt(rho2)),
        window_C2=bool(rho2 <= 4.0),
    )
    P(f"  {tag:<34s} e1={e1:.8g} e2={e2:.8g} rho2={rho2:.10g}  Cmin={d['Cmin']:.6g}"
      f"  mult={d['mult']}  r3={d['r3'] if d['r3'] is None else round(d['r3'],6)}"
      f"  <=4:{'YES' if d['window_C2'] else '**NO**'}  {note}")
    return d


# ---------------------------------------------------------------- 1D chains
P("=" * 100)
P("A1  1D chains  (path = reflecting ends; circle = periodic).  L = 2..6")
P("=" * 100)
res = {}
for L in range(2, 7):
    # path graph Laplacian: 2 on diagonal, -1 off (random-walk generator of the free chain)
    d = np.full(L, 2.0)
    e = np.full(L - 1, -1.0)
    mu = eigh_tridiagonal(d, e, eigvals_only=True)
    res[f"path_L{L}"] = stats(mu, f"path chain L={L}", vacuum=True)
for L in range(2, 7):
    A = np.zeros((L, L))
    for i in range(L):
        A[i, i] = 2.0
        A[i, (i + 1) % L] -= 1.0
        A[i, (i - 1) % L] -= 1.0
    mu = eigh(A, eigvals_only=True)
    res[f"circle_L{L}"] = stats(mu, f"circle chain L={L}", vacuum=True,
                                note="translation invariance -> exact multiplets")
OUT["A1_chains"] = res
flush()

# ---------------------------------------------------------------- box / BC families
P("=" * 100)
P("A2  Particle in a box, ALL four boundary-condition patterns  (L=1)  [T]")
P("    DD/NN/PP tower k^2 -> rho^2=4 ; ND tower ((2k-1)pi/2)^2 -> rho^2=9")
P("=" * 100)
res = {}
for bc in ["DD", "NN", "PP", "ND"]:
    for q in [0.0, 1.0, 10.0]:
        if bc in ("DD", "NN"):
            k = np.arange(1 if bc == "DD" else 1, 5)           # NN drops the k=0 vacuum
            base = (k * np.pi) ** 2
        elif bc == "PP":
            k = np.arange(1, 5)
            base = (2 * k * np.pi) ** 2
        else:
            k = np.arange(1, 5)
            base = ((2 * k - 1) * np.pi / 2) ** 2
        mu = base + q
        res[f"{bc}_q{q:g}"] = stats(mu, f"box {bc} q={q:g}", vacuum=False,
                                    note=f"rho2={mu[1]/mu[0]:.6g}")
OUT["A2_box_BC"] = res
flush()

# ---------------------------------------------------------------- group manifolds
P("=" * 100)
P("A3  One plaquette SU(2) beta=0 = Laplacian on SU(2)~S^3:  mu_l=l(l+2)   [T]")
P("A4  2D YM: mu_R=(g^2/2)C2(R)  (Migdal) -> Casimir tower                 [T]")
P("=" * 100)
res = {}
mu = np.array([l * (l + 2) for l in range(1, 5)], float)
res["plaquette_S3"] = stats(mu, "one plaquette SU(2) (S^3)", vacuum=False,
                            note="mu_l=l(l+2); mult (l+1)^2")
for N in [2, 3, 4, 100]:
    if N == 2:
        mu = np.array([(j * (j + 1)) for j in [0.5, 1.0, 1.5, 2.0]])
    else:
        # SU(N): fund C2=(N^2-1)/(2N); adj C2=N; SU(3) (2,0): C2=10/3
        cf = (N * N - 1) / (2 * N)
        ca = N
        mu = np.array([cf, ca, (10.0 / 3.0 if N == 3 else ca * 1.5)])
    res[f"2dYM_SU{N}"] = stats(mu, f"2D YM SU({N}) Casimir tower", vacuum=False,
                               note="fund -> adj -> next (mult = dim R)")
OUT["A3_plaquette_2dYM"] = res
flush()

# ---------------------------------------------------------------- rotor / 2D U(1)
P("=" * 100)
P("A5  1D rotor and 2D U(1): tower n^2 with EXACT 2-fold (+-n) degeneracy     [T]")
P("    -> raw rho^2 = 1 (Kramers face) ; merged rho^2 = 4")
P("=" * 100)
res = {}
n = np.arange(-6, 7)
mu = (n ** 2).astype(float)
res["rotor_raw"] = stats(mu, "rotor / 2D U(1) (RAW, unmerged +-n)", vacuum=True,
                         merge=False, note="exact 2-fold bottom -> raw rho^2 = 1")
# merged: one n^2 per |n|
mu_m = (np.arange(0, 6) ** 2).astype(float)
res["rotor_merged"] = stats(mu_m, "rotor / 2D U(1) (merged |n|)", vacuum=True)
# numerical: 2D U(1) transfer = heat kernel on U(1), 200-point angle grid
Nt = 200
th = 2 * np.pi * np.arange(Nt) / Nt
dth = th[1] - th[0]
# -d^2/dtheta^2 periodic, 3-point
A = np.zeros((Nt, Nt))
for i in range(Nt):
    A[i, i] = 2 / dth ** 2
    A[i, (i + 1) % Nt] -= 1 / dth ** 2
    A[i, (i - 1) % Nt] -= 1 / dth ** 2
mu_num = eigh(A, eigvals_only=True)
res["rotor_numeric"] = stats(mu_num, "2D U(1) transfer, numeric 200-pt grid", vacuum=True,
                             tol=1e-6, note="lattice check of the n^2 tower")
OUT["A5_rotor_U1"] = res
flush()

# ---------------------------------------------------------------- discs / balls
P("=" * 100)
P("A6  Radial towers (physical Dirichlet Laplacians, q=0, no double well)      [T]")
P("=" * 100)
res = {}
j0 = jn_zeros(0, 3)
j1 = jn_zeros(1, 3)
j2 = jn_zeros(2, 2)
disc_full = np.sort(np.concatenate([j0[:2], j1[:2], j2[:2]]))
res["disc_full"] = stats(disc_full ** 2, "2D disc (full tower)", vacuum=False,
                         note=f"j01={j0[0]:.6f} j11={j1[0]:.6f}")
res["disc_swave"] = stats(j0[:3] ** 2, "2D disc (s-wave/m=0 sector)", vacuum=False,
                          note="same domain, radial sector only -> 5.27 > 4")
# 3D ball: zeros of spherical bessel j_l = J_{l+1/2}
from scipy.optimize import brentq as _bq
from scipy.special import spherical_jn as _sj


def sbez(l, k):
    """first k positive zeros of the spherical Bessel j_l"""
    zs, x = [], 1e-6
    while len(zs) < k:
        x2 = x + 0.2
        if _sj(l, x) * _sj(l, x2) < 0:
            zs.append(_bq(lambda t: _sj(l, t), x, x2))
        x = x2
    return np.array(zs)


sb0 = sbez(0, 3)      # l=0  (= n*pi)
sb1 = sbez(1, 2)      # l=1
P("    spherical Bessel zeros: l=0 ->", np.round(sb0, 6), " l=1 ->", np.round(sb1, 6))
ball_full = np.sort(np.concatenate([sb0[:2], sb1[:2]]))
res["ball_full"] = stats(ball_full ** 2, "3D ball (full tower)", vacuum=False)
res["ball_swave"] = stats(sb0 ** 2, "3D ball (s-wave sector)", vacuum=False,
                          note="spherical s-wave zeros = n*pi -> exactly 4")
OUT["A6_radial"] = res
flush()

# ---------------------------------------------------------------- discrete star
P("=" * 100)
P("A7  Discrete (NON-metric) star K_{1,m}:  spectrum {0, 1^(m-1), m+1}          [T]+[V]")
P("    rho^2 = m+1 > 4 for m>=4  ->  kills C=2 on the combinatorial-graph class")
P("=" * 100)
res = {}
for m in range(2, 7):
    N = m + 1
    A = np.zeros((N, N))
    for j in range(1, N):          # leaf j connected to centre 0
        A[0, 0] += 1
        A[j, j] += 1
        A[0, j] -= 1
        A[j, 0] -= 1
    mu = eigh(A, eigvals_only=True)
    res[f"star_m{m}"] = stats(mu, f"discrete star K_(1,{m})", vacuum=True)
OUT["A7_discrete_star"] = res
flush()

P("=" * 100)
P("A8  summary: models with rho^2 > 4 (C=2 would die on any of these)")
P("=" * 100)
bad = []
for grp in OUT.values():
    for k, v in grp.items():
        if v["window_C2"] is False:
            bad.append((v["tag"], v["rho2"], v["note"]))
for t, r, n in bad:
    P(f"    ** rho^2={r:.6g}  {t}   {n}")
P(f"    total violations: {len(bad)}")
OUT["A8_violations"] = [dict(tag=t, rho2=r, note=n) for t, r, n in bad]
flush()
P("done. json ->", JSON)
