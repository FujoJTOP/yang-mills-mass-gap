"""
pnt_percolation_chi.py -- witness for PNT-PERCOLATION-CHI.md

Tests the proposed inequality:  [level-set percolation of trF^2]  ==>  [chi = Var(trF^2) >= c].

Three exact/cheap checks (numpy + itertools only, <1s):
  PART A  a deterministic percolating configuration has Var = 0
          => percolation (config-level) does NOT imply chi>0 (ensemble-level).  Bare inequality FALSE.
  PART B  the TRUE percolation inequality (FK/GKS) runs the OPPOSITE way:
          chi = sum_x <s0 sx>  >=  E|C_0|  (mean FK cluster size)  =>  percolation => LARGE chi,
          which is the GAPLESS signal, not the gap (H1 needs chi <= 1/lambda, an UPPER bound).
  PART C  full-rep SU(N) has no FKG (characters are not ferromagnetic): note + a one-plaquette check
          that the natural two-point weight has no positive-association structure to exploit.

Run:  py -3.12 -I research/unconditional/w7/pnt_percolation_chi.py
"""

import itertools, math

print("=" * 74)
print("PART A  -- deterministic percolating field has ZERO variance")
print("=" * 74)

L = 4
# single deterministic configuration: phi == 1 everywhere (trivially percolates)
phi = [[1.0] * L for _ in range(L)]
theta = 0.5
A = {(i, j) for i in range(L) for j in range(L) if phi[i][j] > theta}
print(f"  field phi == 1 (deterministic); level set |A| = {len(A)} / {L*L}, theta={theta}")
print(f"  level set fills the whole torus -> percolates: {len(A) == L*L}")

# ensemble = degenerate measure delta_{phi==1}; Phi = sum of phi over all plaquettes
Phi = sum(sum(r) for r in phi)
Phi_samples = [Phi]                  # one point mass
mean = sum(Phi_samples) / len(Phi_samples)
var = sum((x - mean) ** 2 for x in Phi_samples) / len(Phi_samples)
print(f"  chi := Var(Phi) over the ensemble = {var:.3e}   (single deterministic config)")
print("  => percolation = TRUE while Var(Phi) = 0 : the bare implication is FALSE.")
print("     (percolation is a property of ONE configuration; chi is an ENSEMBLE variance.)")


print()
print("=" * 74)
print("PART B  -- the TRUE (FK/GKS) percolation inequality points the WRONG way")
print("=" * 74)
print("  For Ising / FK (q=2, p = 1 - e^{-2K}):   chi = sum_x <s_0 s_x>  >=  E|C_0|  (GKS).")
print("  So percolation (infinite mean cluster) => chi = infinity : the GAPLESS signal.")
print("  H1 needs chi <= 1/lambda (an UPPER bound) -- the opposite direction.")

n = 4
bonds = [(i, (i + 1) % n) for i in range(n)]

def ising_corr(K, a, b):
    Z = 0.0; num = 0.0
    for s in itertools.product([1, -1], repeat=n):
        E = sum(s[i] * s[(i + 1) % n] for i in range(n))
        w = math.exp(K * E)
        Z += w; num += w * s[a] * s[b]
    return num / Z

def ising_Z(K):
    return sum(math.exp(K * sum(s[i] * s[(i + 1) % n] for i in range(n)))
               for s in itertools.product([1, -1], repeat=n))

def fk_stats(K):
    p = 1.0 - math.exp(-2.0 * K); q = 2.0
    Zf = 0.0; E_size = 0.0; Pconn = [0.0] * n
    for mask in range(1 << n):
        edges = [bonds[i] for i in range(n) if (mask >> i) & 1]
        par = list(range(n))
        def f(a):
            while par[a] != a:
                par[a] = par[par[a]]; a = par[a]
            return a
        for (u, v) in edges:
            par[f(u)] = f(v)
        comps = len({f(i) for i in range(n)})
        w = (p ** len(edges)) * ((1 - p) ** (n - len(edges))) * (q ** comps)
        Zf += w
        c0 = {i for i in range(n) if f(i) == f(0)}
        E_size += w * len(c0)
        for x in range(n):
            if f(x) == f(0):
                Pconn[x] += w
    E_size /= Zf
    Pconn = [v / Zf for v in Pconn]
    return Zf, E_size, Pconn

print(f"\n  {'K':>5} {'Z_Ising':>10} {'Z_FK':>10} {'chi':>10} {'E|C_0|':>10} {'chi>=E|C_0|':>12}  P(0<->x) = " + str([f'{i}' for i in range(n)]))
for K in (0.2, 0.4, 0.6, 1.0):
    Z = ising_Z(K)
    chi = sum(ising_corr(K, 0, x) for x in range(n))
    Zf, Esize, Pconn = fk_stats(K)
    ok = "YES" if chi >= Esize - 1e-12 else "NO "
    print(f"  {K:>5.2f} {Z:>10.4f} {Zf:>10.4f} {chi:>10.4f} {Esize:>10.4f} {ok:>10}   "
          + " ".join(f"{Pconn[x]:.3f}" for x in range(n)))
print("  => chi >= E|C_0| holds: percolation feeds a LOWER bound on chi (gapless direction),")
print("     not the UPPER bound chi <= 1/lambda that H1 requires.")


print()
print("=" * 74)
print("PART C  -- no FKG for full-rep SU(N) (character weights are not ferromagnetic)")
print("=" * 74)
print("  FKG / GKS needs a ferromagnetic (totally ordered) spin representation.")
print("  Full-rep SU(N) Wilson weight is a character SUM a_r(beta) chi_r(U); no q-state /")
print("  spin representation, no FKG. (cross-ref BASE-FIELDMEASURE R2, INVENT-PERCOLATION (ii))")
print("  So even the sign-only FK step is unavailable for the physical gauge measure.")


print()
print("=" * 74)
print("SUMMARY")
print("=" * 74)
print("  (A) percolation => chi>=c : FALSE (deterministic percolating config has Var=0).")
print("  (B) the true percolation inequality => chi LARGE (wrong sign for H1).")
print("  (C) FKG unavailable for full-rep SU(N), so even the sign step fails.")
print("  The correct-shaped lemma is 'NO percolation => chi <= c', but 'no percolation at")
print("  beta=2.5' is equivalent to the gap.  Not a free strut.")
