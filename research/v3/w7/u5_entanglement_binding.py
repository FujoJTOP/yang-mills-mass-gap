#!/usr/bin/env python3
"""
U5-ENTANGLEMENT-BINDING witness.  Deterministic, numpy only.
Run:  py -3.12 -I u5_entanglement_binding.py

Question:  can an entanglement / strong-subadditivity / relative-entropy functional give a
rigorous LOWER bound on the atom fraction  w1 = f1^2/M0  (or an UPPER bound on the leakage
1-w1) of the flux-tube state, uniformly through beta=2.5, failing on the free field?

This witness checks the *inequality structure* on solvable models.  It is NOT a proof about
YM.  The four structural facts it illustrates (theorems are in U5-ENTANGLEMENT-BINDING.md):

 (A) SINGLE-MODE EXACT IDENTITY.  For a two-level (single-excitation / single-mode) subsystem
     with reduced DM diag(Z, 1-Z), the atom fraction is w1 = Z and the subsystem entropy is the
     binary entropy S = h(Z) = -Z ln Z - (1-Z) ln(1-Z).  Hence  w1 >= h^{-1}_<(S)  (small root),
     an EXACT lower bound of w1 by a decreasing function of S.  BUT h is SYMMETRIC: h(Z)=h(1-Z),
     so S=0 is compatible with w1=0 (free) AND w1=1 (perfect atom); the bound is VACUOUS
     (w1>=0) at both ends and non-trivial only near Z=1/2.  => ANTICORRELATION: entanglement is
     strong exactly where the atom is small.  This is the acid-test obstruction.

 (B) FREE MASSLESS CHAIN (many-mode Gaussian / free-fermion).  Bipartite entropy S(l) ~ (1/3)log l
     (large, grows), while the free 0++ composite has w1 = 0 exactly.  => LARGE S with ZERO atom:
     a universal floor  w1 >= f(S)  with f(h^{-1}_<) is FALSE (explicit counterexample).

 (C) SSA / MUTUAL INFORMATION.  Mutual information is NEVER negative (SSA).  Random pure states:
     I(A:B)>=0, SSA S(A)+S(C)<=S(AB)+S(BC), Araki-Lieb triangle |S(A)-S(B)|<=S(AB) all hold.
     => "negative mutual information" is impossible; a deficit must be a relative entropy, which
     against a reference tube state is circular (the reference IS the atom).

 (D) BASIS-INDEPENDENCE.  Entanglement S is a functional of rho_A alone.  But w1 is the overlap
     of the source with the OPERATOR-CHOSEN tube direction.  Two operators O,O' on the SAME
     rho_A give different w1 at fixed S:  w1 is not a function of rho_A.  (i.e. no w1 = F(S).)
"""
import numpy as np

np.set_printoptions(precision=6, suppress=True, linewidth=140)
rng = np.random.default_rng(20261124)


def hbin(z):
    """binary entropy, z in [0,1]."""
    z = np.clip(z, 1e-15, 1 - 1e-15)
    return -z * np.log(z) - (1 - z) * np.log(1 - z)


def hbin_inv_left(S):
    """smaller root z* in [0,1/2] of h(z)=S (<= ln2).  h^{-1}_<(S)."""
    ln2 = np.log(2.0)
    S = min(S, ln2)
    lo, hi = 0.0, 0.5
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if hbin(mid) < S:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


print("=" * 82)
print("U5-ENTANGLEMENT-BINDING witness   (numpy %s)" % np.__version__)
print("=" * 82)

# ---------------------------------------------------------------------------
# (A) single-mode exact identity  w1 = Z,  S = h(Z),  w1 >= h^{-1}_<(S)
#     Star model: source site 0 coupled to M bath modes with strength V.
#     ground state weight on the source  = Z.  Single-excitation: rho_0 = diag(Z,1-Z).
# ---------------------------------------------------------------------------
print("\n(A) SINGLE-MODE EXACT IDENTITY   w1 = Z ,   S = h(Z) ,   bound  w1 >= h^{-1}_<(S)")
M = 40                       # bath modes -> Z_min ~ 1/(M+1)
eps_s, eps_b = -1.0, 0.0     # source below bath: V=0 ground = source (Z=1); V large -> Z~1/(M+1)
print(f"    {'V':>8s} {'Z=w1':>9s} {'S=h(Z)':>9s} {'h^{-1}_<(S)':>12s} {'w1>=bound':>10s} {'S(Z)==S(1-Z)':>13s}")
for V in [0.02, 0.10, 0.30, 0.70, 1.20, 2.50, 6.00, 20.0]:
    n = M + 1
    H = np.zeros((n, n))
    H[0, 0] = eps_s
    for i in range(1, n):
        H[i, i] = eps_b
        H[0, i] = H[i, 0] = V          # star coupling source<->bath
    ev, U = np.linalg.eigh(H)
    Z = U[0, 0] ** 2                    # weight of source site in ground state
    S = hbin(Z)
    bnd = hbin_inv_left(S)
    # exact (direct) 2-level reduced entropy of source: -ZlnZ-(1-Z)ln(1-Z)
    assert abs(S - hbin(1 - Z)) < 1e-12, "binary-entropy symmetry must hold"
    print(f"    {V:8.3f} {Z:9.5f} {S:9.5f} {bnd:12.5f} {str(Z>=bnd-1e-12):>10s} {str(abs(S-hbin(1-Z))<1e-12):>13s}")
print("    reading: small S (Z->0 free, Z->1 perfect atom) => bound -> 0 (VACUOUS);")
print("             the bound is strongest only near Z=1/2, the *least* atomic regime.")
print("             => entanglement is ANTICORRELATED with the atom; no floor in the atom regime.")

# other branch: source ABOVE bath -> ground state is the bath, w1 = Z small, S = h(Z) small.
print("    other branch (source above bath -> ground is the bath, w1 = Z small, S small too):")
for eps_s2, V2 in [(5.0, 0.02), (2.0, 0.30), (1.0, 1.0)]:
    n = M + 1
    H = np.zeros((n, n))
    H[0, 0] = eps_s2
    for i in range(1, n):
        H[i, i] = 0.0
        H[0, i] = H[i, 0] = V2
    ev, U = np.linalg.eigh(H)
    Z = U[0, 0] ** 2
    S = hbin(Z)
    print(f"      eps_s={eps_s2:4.1f} V={V2:5.2f}  Z=w1={Z:8.5f}  S={S:8.5f}  bound={hbin_inv_left(S):8.5f}")
print("    => S -> 0 at BOTH w1 -> 0 (free) and w1 -> 1 (perfect atom): S cannot tell them apart.")

# ---------------------------------------------------------------------------
# (B) free massless chain: S(l) ~ (1/3) log l  (LARGE)  while  w1(free composite)=0.
#     correlation-matrix method (free fermions, half filling).
# ---------------------------------------------------------------------------
print("\n(B) FREE MASSLESS CHAIN  ->  S(l) grows ~ (1/3)log l ,  w1(free 0++) = 0 exactly")
N = 401
Hc = np.zeros((N, N))
for i in range(N - 1):
    Hc[i, i + 1] = Hc[i + 1, i] = -1.0
ev, U = np.linalg.eigh(Hc)
occ = U[:, : N // 2]                    # half-filled
C = occ @ occ.T                         # <c_i^dag c_j>
def ent_of_interval(C, lo, hi):
    lam = np.clip(np.linalg.eigvalsh(C[lo:hi, lo:hi]), 1e-12, 1 - 1e-12)
    return -np.sum(lam * np.log(lam) + (1 - lam) * np.log(1 - lam))
print(f"    {'l':>4s} {'S(l)':>9s} {'S/l':>8s} {'(1/3)ln l':>10s}")
mid = N // 2
for l in [8, 16, 32, 64, 128]:
    S = ent_of_interval(C, mid - l // 2, mid + l // 2)
    print(f"    {l:4d} {S:9.4f} {S/l:8.4f} {np.log(l)/3:10.4f}")
print("    reading: S large and GROWING (area-law-with-log), yet the free 0++ composite")
print("             has f1^2 = w1 = 0 (two-particle cut, no atom).  Large S, zero atom.")
print("             => NO universal floor  w1 >= f(S)  with f = h^{-1}_< can exist [counterexample].")

# ---------------------------------------------------------------------------
# (C) SSA / mutual information is never negative  (random pure 4-qubit states)
# ---------------------------------------------------------------------------
print("\n(C) SSA / MUTUAL INFORMATION  (random pure 4-qubit states; 4000 samples)")
def vN(rho):
    w = np.clip(np.linalg.eigvalsh(rho), 1e-15, 1.0)
    return -np.sum(w * np.log(w))
def reduced_from_psi(psi, keep):
    """reduced DM of a 4-qubit pure state `psi` (dim 16) on the sites in `keep`."""
    keep = sorted(keep)
    drop = [s for s in range(4) if s not in keep]
    T = psi.reshape([2] * 4)
    T = np.transpose(T, keep + drop).reshape(2 ** len(keep), 2 ** len(drop))
    return T @ T.conj().T
I_neg = ssa_viol = tri_viol = 0
worst_I = worst_AL = worst_SSA = -1.0
for _ in range(4000):
    psi = rng.standard_normal(16) + 1j * rng.standard_normal(16)
    psi /= np.linalg.norm(psi)
    SA = vN(reduced_from_psi(psi, [0]))
    SB = vN(reduced_from_psi(psi, [1]))
    SAB = vN(reduced_from_psi(psi, [0, 1]))
    SABC = vN(reduced_from_psi(psi, [0, 1, 2]))
    SBCD = vN(reduced_from_psi(psi, [1, 2, 3]))
    SBC = vN(reduced_from_psi(psi, [1, 2]))
    # (1) subadditivity I(A:B)=SA+SB-SAB >= 0
    dI = SA + SB - SAB
    worst_I = max(worst_I, -dI)
    if dI < -1e-9:
        I_neg += 1
    # (2) Araki-Lieb |S(A)-S(B)| <= S(AB)
    dAL = SAB - abs(SA - SB)
    worst_AL = max(worst_AL, -dAL)
    if dAL < -1e-9:
        tri_viol += 1
    # (3) SSA: S(ABC)+S(BCD) >= S(BC)+S(ABCD) ; global pure => S(ABCD)=0
    dSSA = SABC + SBCD - SBC
    worst_SSA = max(worst_SSA, -dSSA)
    if dSSA < -1e-9:
        ssa_viol += 1
print(f"    subadditivity I(A:B)<0 violations : {I_neg}   (min I = {-worst_I:.2e})")
print(f"    Araki-Lieb  |S(A)-S(B)|>S(AB) viol : {tri_viol}   (worst = {worst_AL:.2e})")
print(f"    SSA  S(ABC)+S(BCD)<S(BC) violations: {ssa_viol}   (worst = {worst_SSA:.2e})")
print("    reading: MI>=0 and SSA hold identically. 'Negative mutual information' is impossible;")
print("             an 'entanglement deficit' must be a RELATIVE entropy vs a reference tube,")
print("             which is circular (the reference IS the atom).")

# ---------------------------------------------------------------------------
# (D) basis-independence: w1 is NOT a function of rho_A (same S, different w1)
# ---------------------------------------------------------------------------
print("\n(D) BASIS-INDEPENDENCE   (same rho_sys & S, different operator => different w1)")
a, b, c = 0.10, 0.20, 0.70            # fixed populations on {v, G, perp}
rho_pop = np.array([a, b, c])
S_fixed = -np.sum(rho_pop * np.log(rho_pop))
psi_sys = np.sqrt(rho_pop).astype(complex)          # (sqrt a, sqrt b, sqrt c)
# operator 1: tube direction = |G>  (index 1)
tube1 = np.array([0, 1, 0], complex)
w1_1 = abs(np.vdot(tube1, psi_sys)) ** 2             # = b
# operator 2: tube direction = (|G>+|perp>)/sqrt2
tube2 = np.array([0, 1, 1], complex) / np.sqrt(2)
w1_2 = abs(np.vdot(tube2, psi_sys)) ** 2
print(f"    fixed populations {{v,G,perp}} = {rho_pop}   S_fixed = {S_fixed:.5f}")
print(f"    O_1: tube=|G>              ->  w1 = {w1_1:.5f}")
print(f"    O_2: tube=(|G>+|perp>)/v2  ->  w1 = {w1_2:.5f}")
print(f"    same rho_sys & S,  w1 differs by {abs(w1_1-w1_2):.5f}")
print("    reading: entanglement S is a functional of rho_A ALONE; w1 is the overlap with the")
print("             OPERATOR-CHOSEN tube ray.  w1 is not a function of rho_A => no w1=F(S).")
print("             Any rho_A-only (entanglement) bound is O-blind; the O-dependence is exactly")
print("             the residue f1^2 itself (circular).")

# ---------------------------------------------------------------------------
# (E) 2D-analogue: MI > 0 (correlation / string) while a CONSERVED density has NO atom
#     (the 2D-YM fact trF^2 = H is conserved => <0|trF^2|G> = 0 => w1 = 0, R5 §1.1).
# ---------------------------------------------------------------------------
print("\n(E) 2D-ANALOGUE   correlated state (MI>0) with a CONSERVED density -> w1 = 0")
l = 16
def Sg(Cb):
    lam = np.clip(np.linalg.eigvalsh(Cb), 1e-12, 1 - 1e-12)
    return -np.sum(lam * np.log(lam) + (1 - lam) * np.log(1 - lam))
for gap in [4, 16]:
    A = list(range(mid - l, mid)); B = list(range(mid + gap, mid + gap + l))
    MI = Sg(C[np.ix_(A, A)]) + Sg(C[np.ix_(B, B)]) - Sg(C[np.ix_(A + B, A + B)])
    print(f"    MI(two blocks, gap={gap:2d}) = {MI: .5f}")
print("    conserved density N=sum c^dag c commutes with H => N(t)=N => C(t)=<N^2>-<N>^2 =const")
print("    => KL measure = delta_0 (contact, m=0), NO atom at m>0 => w1 = 0.")
print("    reading: the state is highly correlated (MI>0, like sigma>0 / a confining string),")
print("             yet the trF^2-type conserved density has w1=0 EXACTLY.  Entanglement / MI")
print("             fires where the 0++ atom is zero => wrong observable (1-form/center sector).")

print("\n" + "=" * 82)
print("VERDICT of the witness:  entanglement does NOT give a rigorous binding floor.")
print(" (A) exact single-mode identity but SYMMETRIC -> vacuous where the atom lives (anticorrelation)")
print(" (B) free massless chain: huge S, w1=0 -> no universal floor w1>=f(S)  [counterexample]")
print(" (C) MI>=0 (SSA); 'negative MI' impossible; deficit must be circular relative entropy")
print(" (D) w1 is not a function of rho_A (basis/operator dependence) -> entanglement is O-blind")
print(" The exact missing estimate (see U5-ENTANGLEMENT-BINDING.md) is a COHERENCE / asymmetric")
print(" lower bound, which lands on the polymer residue w1>=w0 = isolation = H*.  Clay 0%.")
print("=" * 82)
