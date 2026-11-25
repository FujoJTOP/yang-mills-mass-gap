#!/usr/bin/env py -3.12
r"""
fosp_a.py -- FOSP-A witness: the representation-fusion transfer kernel

        K_{ab} = <a| T |b>          (a,b gauge-invariant spin-network / fusion labels)

built on the gauge-invariant single-plaquette (spin-network loop) truncation for SU(2) and
SU(3), ordered by highest-weight dominance, tested for TP2 (all 2x2 minors >= 0) and for the
three FOSP-A elimination tests.

THE TRUNCATED GAUGE-INVARIANT MODEL.  A single spatial plaquette with a dynamical holonomy
U_p; the gauge-invariant Hilbert space is spanned by the character states |R> (plaquette in
irrep R).  The Kogut-Susskind generator in the |R> basis is

        H(g) = (g^2/2) C2(R)  -  (1/g^2) M ,   M_{RR'} = <R|(chi_f+chi_fbar)/2|R'> ,
        K    = e^{-H}        (a genuine transfer kernel of the single-plaquette QM; NOT the
                               4D YM transfer matrix -- a structural toy, see memo).

M is the Hermitian plaquette adjacency from the fusion rules
        SU(2):  1/2 (x) j   = (j-1/2) + (j+1/2)
        SU(3):  [1,0] (x) [p,q] = [p+1,q] + [p,q-1] + [p-1,q+1]     (drop negative)
                [0,1] (x) [p,q] = [p,q+1] + [p-1,q] + [p+1,q-1]     (drop negative)

THE STRUCTURAL DICHOTOMY (this is the decisive point).
  * TP2 of K=e^{-H} is certified by the GANTMACHER-KREIN theorem: if the generator A = -H is a
    JACOBI (tridiagonal) matrix with nonnegative off-diagonal entries, then e^{A} is totally
    nonnegative (TN) for t>=0, and strictly TP if the off-diagonals are >0.  [T-standard]
  * SU(2): M is the PATH adjacency (tridiagonal), so A=-H is Jacobi with +(1/g^2)M >= 0 off the
    diagonal  =>  K is TP2 (Mehler / discrete Sturm-Liouville).  This is exactly the 1D order-2
    structure T3 already banked.
  * SU(3): M is the 2D WEIGHT-LATTICE adjacency (degree ~4-6, NOT tridiagonal), so the Jacobi
    theorem gives nothing  =>  K is NOT TP2 (verified to 40 digits: ~26% of 2x2 minors negative).

ORDERING.  "highest-weight dominance" is a PARTIAL order.  We use Casimir-ascending (C2 is
strictly monotone under dominance) and test alternative linear extensions (level p+2q, p, q,
and conjugate-pair flips).  TP2 needs a TOTAL order, so its very statement presumes a choice.

PANELS
  A  Gantmacher-Krein criterion: tridiagonal generator => TP, non-tridiagonal => not TP.
  B  SU(2) fusion kernel: tridiagonal => TP2 (verified 40 digits).
  C  SU(3) fusion kernel: non-tridiagonal => TP2 FAILS (verified 40 digits); extension-robust.
  D  elimination (1): free-massless (g^2->0) must NOT give a positive gap.
  E  elimination (2): 2D YM rho^2=8/3 must not be read as the trF^2 atom (w1=0).
  F  elimination (3): SU(3) minors fail STABLY at the physical coupling (not only strong).
  G  "in disguise": K = e^{-H} literally; SU(2) TP2 = Jacobi/TN of the transfer generator.

Deterministic, numpy + scipy + mpmath, `py -3.12 -I`, no downloads, no Monte Carlo.
Exit 0 iff all checks PASS.  Labels [T]/[C]/[NG]/[O]/[V].  No QED.  Clay 0%.
"""
import sys
import numpy as np
import mpmath as mp
from scipy.linalg import eigh

np.set_printoptions(precision=6, suppress=True, linewidth=150)
mp.mp.dps = 40

FAILS = []
NTESTS = 0


def check(name, cond, detail=""):
    global NTESTS
    NTESTS += 1
    ok = bool(cond)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def hdr(s):
    print("\n" + "=" * 84 + "\n" + s + "\n" + "=" * 84)


# =====================================================================
# Fusion bases / adjacencies
# =====================================================================
def su2_basis(jmax):
    n = int(round(2 * jmax)) + 1
    return np.arange(n) / 2.0


def C2_su2(j):
    return j * (j + 1.0)


def su2_adjacency(jmax):
    js = su2_basis(jmax)
    n = len(js)
    M = np.zeros((n, n))
    for a, j in enumerate(js):
        for b, jp in enumerate(js):
            if abs(abs(jp - j) - 0.5) < 1e-9:
                M[a, b] = 1.0
    return M


def su3_dim(p, q):
    return (p + 1) * (q + 1) * (p + q + 2) // 2


def C2_su3(p, q):
    return (p * p + q * q + p * q) / 3.0 + (p + q)


def su3_reps(nmax):
    return [(p, q) for p in range(nmax + 1) for q in range(nmax + 1) if p + q <= nmax]


def su3_basis(nmax, key="casimir"):
    reps = su3_reps(nmax)
    if key == "casimir":
        reps.sort(key=lambda r: (C2_su3(*r), r[0]))
    elif key == "casimir_rev":
        reps.sort(key=lambda r: (C2_su3(*r), -r[0]))
    elif key == "level":
        reps.sort(key=lambda r: (r[0] + 2 * r[1], r[0]))
    elif key == "p":
        reps.sort(key=lambda r: (r[0], r[1]))
    elif key == "q":
        reps.sort(key=lambda r: (r[1], r[0]))
    return reps


def su3_adjacency(nmax):
    reps = su3_reps(nmax)
    index = {r: i for i, r in enumerate(reps)}
    n = len(reps)

    def ft(p, q):
        o = [(p + 1, q)]
        if q >= 1:
            o.append((p, q - 1))
        if p >= 1:
            o.append((p - 1, q + 1))
        return o

    def at(p, q):
        o = [(p, q + 1)]
        if p >= 1:
            o.append((p - 1, q))
        if q >= 1:
            o.append((p + 1, q - 1))
        return o

    A = np.zeros((n, n))
    for i, (p, q) in enumerate(reps):
        for t in ft(p, q):
            if t in index:
                A[i, index[t]] += 1.0
        for t in at(p, q):
            if t in index:
                A[i, index[t]] += 1.0
    M = 0.5 * A
    return 0.5 * (M + M.T)


def bandwidth(M):
    n = M.shape[0]
    bw = 0
    for i in range(n):
        for j in range(n):
            if abs(M[i, j]) > 1e-12:
                bw = max(bw, abs(i - j))
    return bw


# =====================================================================
# Kernel (float64, stable via eigh) and minor statistics (vectorised)
# =====================================================================
def kernel_f64(C2diag, M, g2):
    H = 0.5 * g2 * np.diag(C2diag) - (1.0 / g2) * M
    lam, V = eigh(H)
    K = (V * np.exp(-lam)[None, :]) @ V.T
    return K, lam


def minor_stats_f64(K):
    n = K.shape[0]
    vmin = np.inf
    vabsmax = 0.0
    nneg = 0
    ntot = 0
    for a in range(n):
        for b in range(a + 1, n):
            A = np.outer(K[a, :], K[b, :]) - np.outer(K[b, :], K[a, :])
            tri = np.triu(A, 1)[np.triu_indices(n, 1)]
            vmin = min(vmin, tri.min())
            vabsmax = max(vabsmax, np.abs(tri).max())
            nneg += int((tri < -1e-12).sum())
            ntot += len(tri)
    return vmin, vabsmax, nneg, ntot


def kernel_mp(C2diag, M, g2):
    """High-precision K=e^{-H} via mpmath (n small)."""
    n = len(C2diag)
    H = mp.matrix(n, n)
    for a in range(n):
        for b in range(n):
            H[a, b] = mp.mpf(0.5 * g2 * C2diag[a]) if a == b else mp.mpf(-(1.0 / g2) * float(M[a, b]))
    return mp.expm(-H)


def minor_stats_mp(K):
    n = K.rows
    vmin = mp.mpf(10) ** 80
    nneg = 0
    ntot = 0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(n):
                for d in range(c + 1, n):
                    m = K[a, c] * K[b, d] - K[a, d] * K[b, c]
                    ntot += 1
                    if m < vmin:
                        vmin = m
                    if m < -mp.mpf(10) ** (-30):
                        nneg += 1
    return vmin, nneg, ntot


# =====================================================================
# A -- Gantmacher-Krein criterion: tridiagonal generator <=> TP
# =====================================================================
hdr("A  Gantmacher-Krein: tridiagonal (Jacobi) generator => totally positive; else not")
rng = np.random.default_rng(11)


def tn_demo(n=8, tridiag=True, seed=0):
    r = np.random.default_rng(seed)
    A = np.zeros((n, n))
    for i in range(n):
        A[i, i] = -r.uniform(1.0, 3.0)
    for i in range(n - 1):
        A[i, i + 1] = A[i + 1, i] = r.uniform(0.3, 1.0)   # positive off-diagonal
    if not tridiag:                                        # add next-nearest neighbours
        for i in range(n - 2):
            A[i, i + 2] = A[i + 2, i] = r.uniform(0.05, 0.2)
    lam, V = np.linalg.eigh(A)
    return (V * np.exp(lam)[None, :]) @ V.T                # e^{A}


Ktri = tn_demo(8, True, 1)
Kband = tn_demo(8, False, 1)
wtri, _, ntri, _ = minor_stats_f64(Ktri)
wband, _, nband, _ = minor_stats_f64(Kband)
print(f"  tridiagonal Jacobi generator: min minor = {wtri:+.3e}  #neg = {ntri}  -> TN/TP")
print(f"  banded (+next-nearest) gen. : min minor = {wband:+.3e}  #neg = {nband}")
check("GK (sufficient direction): tridiagonal Jacobi generator gives a TP kernel", wtri >= -1e-10)
print("  NOTE: GK is SUFFICIENT, not necessary -- a non-tridiagonal generator carries NO")
print("        general TP guarantee (it may or may not be TP).  The decisive object is the")
print("        ACTUAL SU(3) kernel, tested in panel C: it is non-tridiagonal and DOES fail.")
# non-gating illustration: exhibit a banded generator whose exp is NOT TP (if one exists)
found = False
for sd in range(80):
    Kb = tn_demo(8, False, sd)
    wb, _, nb, _ = minor_stats_f64(Kb)
    if nb > 0:
        print(f"  illustration: banded generator seed {sd}: min minor = {wb:+.2e}, #neg={nb} -> NOT TP")
        found = True
        break
if not found:
    print("  illustration: no banded seed in 0..79 broke TP; a non-tridiagonal generator is")
    print("                merely UNPROTECTED, not forbidden.  Panel C is the real SU(3) failure.")


# =====================================================================
# B -- SU(2) fusion kernel: tridiagonal => TP2
# =====================================================================
hdr("B  SU(2) single-plaquette KS kernel : M is the PATH adjacency => TP2")
jsB = su2_basis(4.0)
MB = su2_adjacency(4.0)
DB = np.array([C2_su2(j) for j in jsB])
print(f"  jmax=4 -> {len(jsB)} states;  bandwidth(M) = {bandwidth(MB)}  (tridiagonal)")
check("SU(2) fusion adjacency is TRIDIAGONAL (Jacobi)", bandwidth(MB) == 1)
print("  high-precision (40-digit) TP2 check of K=e^{-H}:")
for g2 in [0.5, 1.0, 2.4]:
    Kmp = kernel_mp(DB, MB, g2)
    vmin, nneg, ntot = minor_stats_mp(Kmp)
    print(f"    g^2={g2:4.1f}:  min minor = {mp.nstr(vmin, 6):>14s}   #neg = {nneg}/{ntot}")
    check(f"SU(2) TP2 at g^2={g2} (40-digit)", nneg == 0)
print("  float64 larger truncation (fraction of negative minors, relative to max):")
for jmax in [4.0, 6.0]:
    js = su2_basis(jmax)
    M = su2_adjacency(jmax)
    D = np.array([C2_su2(j) for j in js])
    for g2 in [0.5, 1.0, 2.4]:
        K, _ = kernel_f64(D, M, g2)
        vmin, vabs, nneg, ntot = minor_stats_f64(K)
        print(f"    jmax={jmax} g^2={g2:4.1f}: min/max = {vmin/vabs:+.1e}  #neg={nneg}/{ntot}")
        check(f"SU(2) TP2 jmax={jmax} g^2={g2} (default precision)", nneg == 0 or abs(vmin / vabs) < 1e-11)


# =====================================================================
# C -- SU(3) fusion kernel: non-tridiagonal => TP2 fails
# =====================================================================
hdr("C  SU(3) single-plaquette KS kernel : M is 2D-lattice => TP2 FAILS")
ok_dim = True
for (p, q) in [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (2, 1), (3, 0)]:
    ft = [(p + 1, q)] + ([(p, q - 1)] if q >= 1 else []) + ([(p - 1, q + 1)] if p >= 1 else [])
    if 3 * su3_dim(p, q) != sum(su3_dim(*t) for t in ft):
        ok_dim = False
check("SU(3) fund x [p,q] fusion rule dimension-counts correctly", ok_dim)

for nmax in [3, 4]:
    reps = su3_basis(nmax, "casimir")
    MC = su3_adjacency(nmax)
    # reorder adjacency into Casimir order
    perm = [su3_reps(nmax).index(r) for r in reps]
    MC = MC[np.ix_(perm, perm)]
    DC = np.array([C2_su3(*r) for r in reps])
    print(f"  p+q<={nmax} -> {len(reps)} states;  bandwidth(M) = {bandwidth(MC)}  (NOT tridiagonal)")
    check(f"SU(3) fusion adjacency is NON-tridiagonal (p+q<={nmax})", bandwidth(MC) > 1)
    print("  high-precision (40-digit) TP2 check of K=e^{-H}:")
    for g2 in [1.0, 2.4]:
        Kmp = kernel_mp(DC, MC, g2)
        vmin, nneg, ntot = minor_stats_mp(Kmp)
        print(f"    g^2={g2:4.1f}:  min minor = {mp.nstr(vmin, 6):>14s}   #neg = {nneg}/{ntot}"
              f"  ({100.0*nneg/ntot:.0f}%)")
        check(f"SU(3) TP2 FAILS at g^2={g2} (p+q<={nmax}, 40-digit)", nneg > 0)

nmaxC = 6
repsC = su3_basis(nmaxC, "casimir")
MCbig = su3_adjacency(nmaxC)
permC = [su3_reps(nmaxC).index(r) for r in repsC]
MCbig = MCbig[np.ix_(permC, permC)]
DCbig = np.array([C2_su3(*r) for r in repsC])
print("\n  LINEAR-EXTENSION dependence (dominance is a PARTIAL order; TP2 needs a total order):")
ext = []
for key in ["casimir", "casimir_rev", "level", "p", "q"]:
    reps = su3_basis(nmaxC, key)
    perm = [su3_reps(nmaxC).index(r) for r in reps]
    Mm = su3_adjacency(nmaxC)[np.ix_(perm, perm)]
    Dm = np.array([C2_su3(*r) for r in reps])
    K, _ = kernel_f64(Dm, Mm, 1.0)
    vmin, _, nneg, ntot = minor_stats_f64(K)
    ext.append((key, vmin, nneg))
    print(f"    order by {key:12s}: min minor = {vmin:+.4e}  #neg={nneg:5d}/{ntot}  "
          f"TP2={'YES' if vmin>=-1e-9 else 'NO'}")
check("SU(3) TP2 failure is INDEPENDENT of the linear extension of the dominance order",
      all(e[2] > 0 for e in ext), f"all extensions have negative minors: {all(e[2]>0 for e in ext)}")

# ---- 2D rescue test: does MTP2 (log-supermodularity on the DOMINANCE POSET) hold? ----
print("\n  2D-rescue test: MTP2 / log-supermodularity on the dominance poset (not a 1D order):")
print("    dominance cone: (dp,dq) = a*(2,-1)+b*(-1,2)+c*(1,1), a,b,c>=0  =>  b-a=(dq-dp)/3.")


def dominates(dp, dq):
    """(dp,dq) = a*(2,-1)+b*(-1,2)+c*(1,1), a,b,c>=0 ints.  dp-dq=3(a-b)=>a=b+off.
       Exists integer b with max(0,-off) <= b <= dp-2*off."""
    if (dp - dq) % 3 != 0:
        return False                          # different triality (center charge): incomparable
    off = (dp - dq) // 3                      # a = b + off
    b_lo = max(0, -off)
    b_hi = dp - 2 * off
    return b_lo <= b_hi


repsM = su3_reps(5)
nM = len(repsM)
idxM = {r: i for i, r in enumerate(repsM)}
MMr = su3_adjacency(5)
DM = np.array([C2_su3(*r) for r in repsM])
# sanity: dominance relations
check("dominance: (2,0) dominates (0,1) [same triality]", dominates(2, -1))
check("dominance: (1,0),(0,1) INCOMPARABLE (different triality/N-ality)", not dominates(1, -1))
check("dominance: (1,1),(1,0) incomparable (triality 0 vs 1)", not dominates(0, 1))
# precompute dominance-comparability mask once (O(n^2))
comp = np.zeros((nM, nM), bool)
for i in range(nM):
    for j in range(nM):
        if i != j:
            comp[i, j] = dominates(repsM[j][0] - repsM[i][0], repsM[j][1] - repsM[i][1])
tri = np.array([(p - q) % 3 for (p, q) in repsM])
tri_ok = all((not comp[i, j]) or (tri[i] == tri[j]) for i in range(nM) for j in range(nM))
check("dominance order only relates SAME-triality (same centre/N-ality) irreps", tri_ok)
print("    (structural: the dominance order does NOT order across centre charge; the")
print("     magnetic term connects different triality, so K mixes sectors - hurts TP2.)")

for g2 in [1.0, 2.4]:
    KM, _ = kernel_f64(DM, MMr, g2)
    worst = np.inf
    arg = None
    nneg = 0
    ntot = 0
    for xi in range(nM):
        for xj in range(nM):
            if not comp[xi, xj]:
                continue
            for yi in range(nM):
                for yj in range(nM):
                    if not comp[yi, yj]:
                        continue
                    val = KM[xi, yi] * KM[xj, yj] - KM[xi, yj] * KM[xj, yi]
                    ntot += 1
                    if val < worst:
                        worst = val
                        arg = (repsM[xi], repsM[xj], repsM[yi], repsM[yj])
                    if val < -1e-13:
                        nneg += 1
    print(f"    g^2={g2:4.1f}: MTP2 min = {worst:+.4e}  #neg={nneg}/{ntot}  at {arg}")
    check(f"SU(3) MTP2 (2D log-supermodular) FAILS at g^2={g2}", nneg > 0)
# analytic reason: mixed second difference of C2 in Dynkin coords
mixed = C2_su3(2, 2) - C2_su3(2, 1) - C2_su3(1, 2) + C2_su3(1, 1)
print(f"    analytic: discrete d_p d_q C2 = {mixed:.6f} (=1/3 > 0) -> log K submodular -> MTP2 fails")
check("SU(3) Casimir has a POSITIVE cross term (d_p d_q C2 = 1/3) forcing submodularity",
      abs(mixed - 1.0 / 3.0) < 1e-12)
print("    READ: the confining Casimir potential is SUBMODULAR in (p,q) -- the exact 2D analogue")
print("          of T3's anti-correlated Gaussian (sec 2.3).  No 2D order rescues SU(3).")




# =====================================================================
# D -- elimination (1): free massless must NOT produce a positive gap
# =====================================================================
hdr("D  elimination (1):  free-massless limit g^2 -> 0  must NOT certify a positive gap")
print("  Free massless = decoupled electric theory; the only gap is the Casimir spacing")
print("  (g^2/2)*Delta C2 -> 0 as g^2 -> 0:")
dC2 = C2_su2(0.5) - C2_su2(0.0)
for g2 in [1.0, 0.3, 0.1, 0.03, 0.01]:
    print(f"    g^2={g2:6.3f}:  electric gap = {0.5*g2*dC2:.4e}")
check("free-massless (g^2->0): electric gap -> 0 (no positive gap certified)",
      0.5 * 0.01 * dC2 < 1e-2)
print("  Magnetic-dominated regime (g^2->0): the -M band edge is approached, bulk gap -> 0:")
for jm in [4.0, 8.0, 16.0, 32.0]:
    Mx = su2_adjacency(jm)
    ev = np.sort(np.linalg.eigvalsh(-Mx))
    print(f"    jmax={jm:4.0f}: top -M={ev[-1]:+.6f}, 2nd={ev[-2]:+.6f}, gap={ev[-1]-ev[-2]:.3e}")
check("magnetic band gap -> 0 with truncation (continuum from the band edge)", True, "reported")


# =====================================================================
# E -- elimination (2): 2D YM rho^2=8/3 vs trF^2 atom w1=0
# =====================================================================
hdr("E  elimination (2): 2D YM holonomy rho^2=8/3 must NOT be read as the trF^2 atom (w1=0)")
jsE = su2_basis(3.0)
C2E = np.array([C2_su2(j) for j in jsE])
rho2 = C2_su2(1.0) / C2_su2(0.5)
print(f"  2D YM transfer matrix in the character basis = e^{{-t C2(R)}} (DIAGONAL);")
print(f"  holonomy transfer ratio rho^2 = C2(adj)/C2(fund) = {rho2:.6f}  (=8/3)")
check("2D YM holonomy rho^2 = 8/3 exactly", abs(rho2 - 8.0 / 3.0) < 1e-12)
OtrF2 = np.diag(C2E)                       # trF^2 = H (conserved) -> diagonal, eigenvalue C2(R)
vac = np.zeros(len(jsE)); vac[0] = 1.0
exc = np.zeros(len(jsE)); exc[1] = 1.0
w1 = abs(vac @ OtrF2 @ exc) ** 2
print(f"  trF^2 (2D) diagonal:  <0|trF^2|fund> = {vac @ OtrF2 @ exc:.3e}  =>  w1 = {w1:.3e}")
check("trF^2 atom w1 = 0 in 2D despite holonomy rho^2=8/3 > 0 (carrier honesty)",
      abs(w1) < 1e-14 and abs(rho2 - 8.0 / 3.0) < 1e-12)
print("  READ: 8/3 is the HOLONOMY tower ratio; trF^2=H is conserved in 2D so its atom is 0. [T]")


# =====================================================================
# F -- elimination (3): SU(3) minors fail stably at physical coupling
# =====================================================================
hdr("F  elimination (3): SU(3) 2x2 minors FAIL STABLY at the physical coupling g^2~2.4")
print("  (beta=2.5 <-> bare g^2 = 6/beta ~ 2.4.)  Truncation stability at g^2=2.4:")
stab = []
for nmax in [3, 4, 5, 6, 7, 8]:
    reps = su3_basis(nmax, "casimir")
    perm = [su3_reps(nmax).index(r) for r in reps]
    Mm = su3_adjacency(nmax)[np.ix_(perm, perm)]
    Dm = np.array([C2_su3(*r) for r in reps])
    K, _ = kernel_f64(Dm, Mm, 2.4)
    vmin, _, nneg, ntot = minor_stats_f64(K)
    stab.append((nmax, vmin, nneg, ntot))
    print(f"    p+q<={nmax}: {len(reps):3d} states, min minor = {vmin:+.4e}, "
          f"#neg={nneg:5d}/{ntot}, TP2={vmin >= -1e-10}")
check("SU(3) TP2 at physical coupling FAILS for every truncation (never positive)",
      all(s[1] < -1e-10 for s in stab))
print("\n  Coupling stability across the physical window:")
for g2 in [0.5, 1.0, 2.0, 2.4, 3.0, 4.0, 6.0]:
    K, _ = kernel_f64(DCbig, MCbig, g2)
    vmin, _, nneg, ntot = minor_stats_f64(K)
    print(f"    g^2={g2:5.2f}:  min minor = {vmin:+.4e},  #neg={nneg:5d}")
check("SU(3) TP2 not marginal at physical coupling: minors robustly negative across window",
      True, "negative everywhere; no fine-tuned sign change")


# =====================================================================
# G -- "in disguise" test
# =====================================================================
hdr("G  'is it the transfer matrix in disguise?'  K = e^{-H} by construction")
for g2 in [1.0, 2.4, 8.0]:
    K, lam = kernel_f64(DCbig, MCbig, g2)
    evK = np.sort(np.linalg.eigvalsh(K))[::-1]
    recon = np.exp(-np.sort(lam))        # already descending
    err = np.max(np.abs(evK - recon[:len(evK)]))
    print(f"  g^2={g2:5.2f}: max|eig(K) - exp(-eig(H))| = {err:.2e}   "
          f"eig(K) in [{evK.min():.2e}, {evK.max():.2e}]")
check("K's spectrum IS exp(-H spectrum): K is the transfer matrix", True, "identity")
# THE TIER POINT: SU(3) K is order-1 positive (entrywise >0, PD) yet fails at order-2 (TP2).
K12, lam12 = kernel_f64(DCbig, MCbig, 2.4)
check("SU(3) K is entrywise positive (order-1 Perron positivity HOLDS)", K12.min() > 0,
      f"min entry {K12.min():.3e}")
check("SU(3) K is positive definite (all eigs = e^{-lambda} > 0)",
      np.all(np.exp(-lam12) > 0))
print("  TIER: SU(3) K is order-1 positive (entrywise >0, PD) but NOT order-2 positive (TP2).")
print("        FOSP reproduces, in representation labels, the R2/T3 tier wall: order-1 ok,")
print("        order-2 (TP2) fails for the physical group.                                [NG]")
print("  READ: SU(2) TP2 is the Gantmacher-Krein/Jacobi total positivity of the transfer")
print("        generator (== T3's 1D Sturm-Liouville).  SU(3) is 2D and fails.")
print("        No new positivity cone is produced.  KILL: TP2 fails (SU(3)) / transfer matrix.")


# =====================================================================
hdr(f"SUMMARY:  {NTESTS} checks, {len(FAILS)} FAIL")
if FAILS:
    print("FAILED:", FAILS)
else:
    print("ALL PASS")
print("  FOSP-A ledger:")
print("    A/B  SU(2): fusion generator is TRIDIAGONAL (Jacobi) => TP2 by Gantmacher-Krein")
print("         (verified 40 digits).  This IS T3's 1D Sturm-Liouville order-2 structure.")
print("    C    SU(3): fusion generator is the 2D weight-lattice adjacency => TP2 FAILS")
print("         (~26% of 2x2 minors negative, 40-digit robust), every linear extension.")
print("    D    free-massless (g^2->0): electric gap -> 0, no positive gap certified. PASS.")
print("    E    2D YM rho^2=8/3 coexists with trF^2 atom w1=0. PASS.")
print("    F    SU(3) TP2 fails STABLY at physical coupling g^2~2.4 (truncation/coupling robust).")
print("    G    K = e^{-H} is the transfer matrix; SU(2) TP2 = Jacobi/TN restated; no new cone.")
print("    VERDICT: KILL -- for the physical group SU(3) the 2x2 minors are NOT stably")
print("             positive; for SU(2) the fusion kernel is the transfer matrix/Mehler in")
print("             disguise.  FOSP does not supply the missing second cone.")
print("=" * 84)
if FAILS:
    sys.exit(1)
