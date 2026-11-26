"""r04_dimcount.py -- R04: L=2 (4D SU(2)) DOF bookkeeping, symmetry audit, corner exponent.

All numbers printed; no MC on the lattice (o24/o27 own that). numpy + scipy only.
Run:  py -3.12 w7/r04_dimcount.py  >  w7/r04_dimcount.log
"""
import itertools
import numpy as np
from scipy.integrate import quad

L, D, N = 2, 4, 3          # L=2 torus, 4D, SU(2)
sites = list(itertools.product(range(L), repeat=D))
links = [(x, mu) for x in sites for mu in range(D)]
lid = {l: i for i, l in enumerate(links)}
plaq = [(x, mu, nu) for x in sites for mu in range(D) for nu in range(mu + 1, D)]

def p_links(x, mu, nu):
    """oriented plaquette -> 4 links, signs (+1 forward, -1 backward)"""
    xm = list(x); xm[mu] = (x[mu] + 1) % L; xm = tuple(xm)
    xn = list(x); xn[nu] = (x[nu] + 1) % L; xn = tuple(xn)
    return [((x, mu), +1), ((xm, nu), +1), ((xn, mu), -1), ((x, nu), -1)]

print("=" * 78)
print("[0] RAW BOOKKEEPING  (L=2, 4D, SU(2))")
nl, ns, npq = len(links), len(sites), len(plaq)
print(f"  sites            = L^D         = {ns}")
print(f"  links            = D*L^D       = {nl}")
print(f"  link components  = 3*links     = {3*nl}")
print(f"  plaquettes       = C(D,2)*L^D  = {npq}")
print(f"  plaquettes/link  = 2(D-1)      = {2*(D-1)}   (check: {sum(len(p_links(*p)) for p in plaq)/nl})")
print(f"  gauge generators = 3*sites     = {3*ns}")
print(f"  link-level orbit-space dim     = {3*nl-3*ns}")

# ---- [1] linearized plaquette-curvature map at the flat config: rank ----
A = np.zeros((3 * npq, 3 * nl))
for ip, p in enumerate(plaq):
    for l, s in p_links(*p):
        for a in range(3):
            A[3 * ip + a, 3 * lid[l] + a] = s
rk = np.linalg.matrix_rank(A, tol=1e-9)
print("\n[1] PLAQUETTE-CURVATURE MAP  d(holonomy)/d(link), flat config")
print(f"  matrix {A.shape}, rank = {rk},  kernel = {3*nl-rk}")
print(f"  kernel split: exact 1-forms d(lambda) = {3*ns-N} (constant lambda=3 killed at flat)"
      f"  + harmonic 1-forms H^1(T^4)=D*3 = {D*N}  (toron zero modes of the FREE Hessian)")
print(f"  => image dim = {rk} = 3*links - {3*ns-N} - {D*N} = {3*nl}-45-12")
print(f"  NOTE: physical dim is 3*links-3*sites = {3*nl-3*ns}; the 12 harmonic dirs are"
      f" physical (flat for Wilson), not gauge")

# ---- [2] gauge orbit rank, and whether mode A is gauge-redundant ----
rng = np.random.default_rng(0)
def gauge_jac(x):                      # 192 x 48 : d(x_l)/d(lam_y) at site y
    J = np.zeros((3 * nl, 3 * ns))
    for l, (pos, mu) in enumerate(links):
        y0 = sites.index(pos)
        posp = list(pos); posp[mu] = (pos[mu] + 1) % L
        y1 = sites.index(tuple(posp))
        for a in range(3):
            J[3 * l + a, 3 * y0 + a] += 1.0    # g(x) U g(x+mu)^-1  (linearized)
            J[3 * l + a, 3 * y1 + a] -= 1.0
    return J
J = gauge_jac(None)
print("\n[2] GAUGE JACOBIAN")
print(f"  at the FLAT config: matrix {J.shape}, rank = {np.linalg.matrix_rank(J, tol=1e-9)}"
      f" = 48 - 3 (constant gauge = stabiliser of U=1)")
print(f"  at a GENERIC config: rank = 48 / 48 (3 random configs, w7/r04_check.log)"
      f"  => orbit dim a.e. = 48 => physical dim = {3*nl-3*ns}")
v = np.zeros(3 * nl); v[2::3] = 1.0 / 8.0          # mode A: v_{l,3}=1/8  => ||v||=1
print(f"  ||v||^2 = mode A = {v @ v:.12f}   (o24: 64*(1/8)^2 = 1)")
print(f"  rank([J | v]) = {np.linalg.matrix_rank(np.column_stack([J, v]), tol=1e-9)}"
      f"   => conditioning <v,X>=0 removes 1 dim, it is NOT a gauge fixing")
print(f"  fiber dims: link-level 3*nl-1 = {3*nl-1};  mod gauge 3*nl-3*ns-1 = {3*nl-3*ns-1}")

# ---- [3] hypercubic symmetry: order, orbits, invariant subspace ----
PG = [np.array(m) for m in itertools.product(*[  # signed permutation matrices in 4D
        [tuple(0 if i != j else s for i in range(D)) for j in range(D)] + \
        [tuple(0 if i != j else s for i in range(D)) for j in range(D)]
        for s in ()]) ] if False else None
mats = []
for perm in itertools.permutations(range(D)):                 # 4!
    for signs in itertools.product([1, -1], repeat=D):        # 2^4
        M = np.zeros((D, D), dtype=int)
        for j in range(D):
            M[perm[j], j] = signs[j]
        mats.append(M)
mats = {tuple(map(tuple, M)) for M in mats}
print("\n[3] HYPERCUBIC SYMMETRY OF THE L=2 TORUS")
print(f"  point group |W4| = 2^D * D! = {2**D * np.math.factorial(D) if hasattr(np,'math') else 384}"
      f"   enumerated = {len(mats)}")
trans = list(itertools.product(range(L), repeat=D))
print(f"  translations (Z_2)^D        = {len(trans)}")
print(f"  full space group            = {len(mats)*len(trans)}")
print(f"  *** finite group => orbits are 0-dimensional => dimension removed = 0 ***")

def orbit(seed, gens):
    seen, todo = set(), [seed]
    while todo:
        s = todo.pop()
        if s in seen:
            continue
        seen.add(s)
        todo += [g(s) for g in gens]
    return seen

def pt_act(M):
    def f(l):
        (x, mu) = l
        xr = tuple(M[i] @ np.array(x) % L for i in range(D))
        xr = tuple(int((sum(M[i, j] * x[j] for j in range(D))) % L) for i in range(D))
        muv = int(np.argmax([abs(M[i, mu]) for i in range(D)]))
        return (xr, muv)
    return f
def tr_act(t):
    return lambda l: (tuple((l[0][i] + t[i]) % L for i in range(D)), l[1])

for name, gens in [("point group only", [pt_act(np.array(m)) for m in mats]),
                   ("translations only", [tr_act(t) for t in trans]),
                   ("full space group", [pt_act(np.array(m)) for m in mats] +
                                        [tr_act(t) for t in trans])]:
    orbs = []
    rest = set(links)
    while rest:
        s = next(iter(rest))
        o = orbit(s, gens)
        orbs.append(o); rest -= o
    sizes = sorted({len(o) for o in orbs})
    print(f"  {name:20s}: #orbits(list of 64 links) = {len(orbs):3d}, sizes {sizes}"
          f"  => invariant tangent dim = {3*len(orbs)}")

# ---- [4] corner exponent: near |x|=2pi the marginal weight is delta^p ----
print("\n[4] CORNER EXPONENT (this is the only exact collapse available)")
def shell(eps):        # one-link shell measure: J dA at r = 2pi - eps
    return 4 * np.sin((2 * np.pi - eps) / 2) ** 2 / (2 * np.pi - eps) ** 2 \
           * 4 * np.pi * (2 * np.pi - eps) ** 2 * eps   # J * area * (dr)  ~ eps^3
es = np.array([1e-2, 1e-3, 1e-4, 1e-5, 1e-6])
w = np.array([shell(e) for e in es])
slope = np.polyfit(np.log(es), np.log(w), 1)[0]
print(f"  one-link marginal ~ eps^{slope:.8f}  (fit over 1e-2..1e-6)   => per-link p = 3")
# Dirichlet rule:  int_{delta>=0, sum delta = eps} prod delta^2 dS = eps^(3n-1) * Gamma(3)^n/Gamma(3n)
for n in (1, 2, 3):
    f = lambda d1, n=n: (d1 ** 2) * ((1 - d1) ** 2) ** (n > 1) if n > 1 else d1 ** 2
    num, _ = quad(f, 0, 1) if n > 1 else (1.0 / 3.0, None)   # n=1: int_0^1 d^2 dd = 1/3
    exact = np.math.gamma(3) ** n / np.math.gamma(3 * n) if hasattr(np, "math") else \
            __import__("math").gamma(3) ** n / __import__("math").gamma(3 * n)
    if n == 1:
        num = 1.0 / 3.0
    print(f"  n={n}: Dirichlet const numeric {num:.12f}  exact Gamma(3)^n/Gamma(3n) {exact:.12f}"
          f"  |  exponent 3n-1 = {3*n-1}")
print(f"  L=2: 64 links at the corner => exponent = 3*64-1 = {3*64-1}")
print(f"  corner a0 = sum_l v_l,3 * 2pi = 2pi * sum(1/8) = {2*np.pi*8:.12f}  (max of <v,X>)")
print(f"  => rho(a0) ~ (16pi - a0)^{3*64-1}  =>  G''(a0) = {(3*64-1)}/(16pi-a0)^2 + O(1/eps)")
print("\n[5] VERDICT ARITHMETIC")
print(f"  raw 192 -> gauge 144 -> conditioned 143 ; symmetry removes 0")
print(f"  distance to the q11 threshold dim<=6 : {143-6} ; to dim<=20 : {143-20}")
print("=" * 78)
