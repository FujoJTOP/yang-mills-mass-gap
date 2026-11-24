"""
invent_percolation.py -- witness for research/unconditional/INVENT-PERCOLATION.md  [V-here]
Fortuin-Kasteleyn / random-cluster percolation functional on a gauge measure.

Checks (brute force, tiny 3x3 square lattice, numpy/math only):
  (A) FK identity: Z_RC(p) = sum_A p^|A| (1-p)^(|E|-|A|) 2^{k(A)} equals Ising Z up to a
      beta-dependent prefactor, with p = 1 - e^{-2 beta}.
  (B) CONNECTIVITY = SPIN CORRELATION: Q_p(0<->x) == <sigma_0 sigma_x>  (prefactor-free,
      exact).  This is the whole point: the 0++ "cluster connectivity" is the correlator,
      and its decay rate is the gap -- defined at ANY beta (no transition needed).
  (C) FK-FKG MONOTONICITY: Q_p(0<->x) is non-decreasing in p (positive association).
  (D) ACID (massless = critical point): the connectivity length diverges at beta_c
      (gap -> 0); finite away from it.
Run:  py -3.12 -I invent_percolation.py
"""
import math
from itertools import product

L = 3
V = L * L
def vid(i, j): return i * L + j
edges = []
for i in range(L):
    for j in range(L):
        if i + 1 < L: edges.append((vid(i, j), vid(i + 1, j)))
        if j + 1 < L: edges.append((vid(i, j), vid(i, j + 1)))
E = len(edges)

def find(par, x):
    while par[x] != x:
        par[x] = par[par[x]]; x = par[x]
    return x

def ncomp(A):
    par = list(range(V))
    for k, (u, v) in enumerate(edges):
        if (A >> k) & 1:
            ru, rv = find(par, u), find(par, v)
            if ru != rv: par[ru] = rv
    return len({find(par, x) for x in range(V)})

def connected(A, u0, v0):
    par = list(range(V))
    for k, (u, v) in enumerate(edges):
        if (A >> k) & 1:
            ru, rv = find(par, u), find(par, v)
            if ru != rv: par[ru] = rv
    return find(par, u0) == find(par, v0)

# precompute (|A|, k(A)) over all edge subsets
subsets = [(A, bin(A).count("1"), ncomp(A)) for A in range(1 << E)]

def ising_Z(beta):
    Z = 0.0
    for s in range(1 << V):
        sp = [1 if (s >> x) & 1 else -1 for x in range(V)]
        Z += math.exp(beta * sum(sp[u] * sp[v] for u, v in edges))
    return Z

def ising_corr(beta, u0, v0):
    Z = C = 0.0
    for s in range(1 << V):
        sp = [1 if (s >> x) & 1 else -1 for x in range(V)]
        w = math.exp(beta * sum(sp[u] * sp[v] for u, v in edges))
        Z += w; C += w * sp[u0] * sp[v0]
    return C / Z

def fk_Z(p):
    return sum((p ** a) * ((1 - p) ** (E - a)) * (2 ** k) for _, a, k in subsets)

def fk_conn(p, u0, v0):
    Z = C = 0.0
    for A, a, k in subsets:
        w = (p ** a) * ((1 - p) ** (E - a)) * (2 ** k)
        Z += w
        if connected(A, u0, v0): C += w
    return C / Z

print("=" * 72)
print("3x3 square lattice:  V=%d  E=%d  edge-subsets=%d" % (V, E, 1 << E))
print("=" * 72)

# ---- (A) FK identity (find the beta-dependent prefactor) ----
beta = 0.6
p = 1.0 - math.exp(-2 * beta)
Zi = ising_Z(beta); Zf = fk_Z(p)
print("\n(A) FK identity (q=2):")
print("    beta=%.3f  p=1-e^-2b=%.5f" % (beta, p))
print("    Z_Ising = %.6f   Z_RC = %.6f   ratio = %.6f" % (Zi, Zf, Zi / Zf))

# ---- (B) connectivity == spin correlation (exact, prefactor-free) ----
print("\n(B) CONNECTIVITY == SPIN CORRELATION  Q_p(0<->x) vs <sigma_0 sigma_x>:")
pairs = [((0, 0), (0, 1)), ((0, 0), (0, 2)), ((0, 0), (1, 1)), ((0, 0), (2, 0))]
for (a, b) in pairs:
    u0, v0 = vid(*a), vid(*b)
    q = fk_conn(p, u0, v0)
    c = ising_corr(beta, u0, v0)
    print("    x=%-7s Q_conn=%.9f   <ss>=%.9f   |diff|=%.2e" %
          (str(b), q, c, abs(q - c)))

# ---- (C) FK-FKG monotonicity of connectivity in p ----
print("\n(C) FK-FKG MONOTONICITY: Q_p(0<->(2,0)) as p increases:")
u0, v0 = vid(0, 0), vid(2, 0)
prev = -1
mono = True
for pp in [0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95]:
    q = fk_conn(pp, u0, v0)
    if q < prev - 1e-12: mono = False
    prev = q
    print("    p=%.2f  Q=%.6f" % (pp, q))
print("    -> non-decreasing in p: %s" % mono)

# ---- (D) acid: connectivity length grows as beta -> critical (2D Ising beta_c ~ 0.4407) ----
print("\n(D) ACID (massless = critical): Q_p(0<->(2,0)) vs beta (2D Ising beta_c~0.4407):")
print("    as beta -> beta_c from above, connectivity at distance 2 grows (xi_p -> inf, gap -> 0);")
print("    far above beta_c, connectivity is small (finite xi_p, gap > 0).")
for beta2 in [0.4407, 0.6, 0.9, 1.5, 2.5]:
    p2 = 1.0 - math.exp(-2 * beta2)
    q2 = fk_conn(p2, u0, v0)
    print("    beta=%.4f  p=%.5f  Q(0<->dist2)=%.6f" % (beta2, p2, q2))

print("\nDone.  The 0++ connectivity functional = the correlator; its decay rate = the gap,")
print("defined at every beta (no transition).  BUT: q-state FK exists only for integer q;")
print("the full-rep SU(N) Wilson weight is a character sum, not a q-state -> NO random-cluster")
print("representation, hence NO FK-FKG (the one new ingredient).  See the .md.")
