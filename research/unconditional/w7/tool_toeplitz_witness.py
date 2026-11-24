"""
tool_toeplitz_witness.py -- TOOL-TOEPLITZ witness.

Question: can a Christoffel / Szego (Toeplitz-Hankel) construction give a NON-INTEGRAL
LOWER bound on the KL atom weight f1^2 that escapes the OPE-RHO moment-blindness?

We test the two standard pieces:
  (A) Christoffel function  lambda_n(x) = 1/(v^T H^{-1} v),  H_ij = M_{i+j},  v=(1,x,..,x^{n-1})
      -> classical UPPER bound on the atom:  f1^2 <= lambda_n(x),  lambda_n DOWN to f1^2.
  (B) Gauss / Stieltjes  w_1^{(n)} (smallest-node weight of the n-node Gauss rule)
      -> the would-be LOWER bound:  w_1^{(n)} UP to f1^2 (direction tested).

Both are functions of finitely many moments M_0..M_{2n-2} (Hankel), i.e. INTEGRAL data:
a continuum matching those moments gives the SAME value -> finite-n blind (OPE-RHO).

Run: py -3.12 -I research/unconditional/w7/tool_toeplitz_witness.py
"""
import numpy as np
np.set_printoptions(precision=6, suppress=True, linewidth=120)


def moments(xs, ws, K):
    return np.array([np.sum(ws * xs ** k) for k in range(K)])


def christoffel(xs, ws, x, n):
    """lambda_n(x) = 1/(v^T H^{-1} v); uses ONLY M_0..M_{2n-2}."""
    M = moments(xs, ws, 2 * n - 1)            # M_0 .. M_{2n-2}
    H = np.array([[M[i + j] for j in range(n)] for i in range(n)])
    v = np.array([x ** k for k in range(n)])
    Kn = float(v @ np.linalg.solve(H, v))
    return 1.0 / Kn


def gauss(xs, ws, n):
    """n-node Gauss rule from the measure (orthonormal polys + Jacobi). Returns nodes, weights.
    The rule is exact for deg <= 2n-1, i.e. determined by M_0..M_{2n-1} (integral data)."""
    # orthonormal polynomials via modified Gram-Schmidt with double reorth
    p = []
    for k in range(n):
        v = xs.astype(float) ** k
        for _ in range(2):
            for q in p:
                v = v - np.sum(ws * v * q) * q
        p.append(v / np.sqrt(np.sum(ws * v * v)))
    alpha = np.array([np.sum(ws * xs * p[k] * p[k]) for k in range(n)])
    beta = np.zeros(n)
    for k in range(1, n):
        beta[k] = (np.sum(ws * xs * p[k - 1] * p[k])) ** 2
    J = np.diag(alpha)
    for k in range(1, n):
        J[k, k - 1] = J[k - 1, k] = np.sqrt(beta[k])
    evals, evecs = np.linalg.eigh(J)
    wts = np.sum(ws) * evecs[0, :] ** 2
    return evals, wts


def make_atom_continuum(atom_w, atom_x, cont_w, lo, hi, ncont=400):
    xc = np.linspace(lo, hi, ncont)
    xs = np.concatenate(([atom_x], xc))
    ws = np.concatenate(([atom_w], np.full(ncont, cont_w / ncont)))
    return xs, ws


print("=" * 74)
print("MEASURE A: atom f1^2=0.5 at s=1, + continuum 0.5 on [4,6]  (gap 3)")
xa, wa = make_atom_continuum(0.5, 1.0, 0.5, 4.0, 6.0)

print("\n(A) Christoffel lambda_n(1)  [expect DOWN to 0.5, an UPPER bound]")
for n in range(1, 11):
    print(f"   n={n:2d}: lambda_n(1) = {christoffel(xa, wa, 1.0, n):.6f}")

print("\n(B) Gauss smallest-node weight w_1^(n)  [direction?]  and smallest node")
for n in range(1, 11):
    nodes, wts = gauss(xa, wa, n)
    i = np.argmin(nodes)
    print(f"   n={n:2d}: node_1 = {nodes[i]:.5f}   w_1 = {wts[i]:.6f}   (sum w = {wts.sum():.5f})")

print("\n(C) blindness: lambda_n uses only M_0..M_{2n-2}. Two DIFFERENT continua with the")
print("    SAME first 4 moments must give the SAME lambda_3(1) and the same 2-node Gauss.")
xb, wb = make_atom_continuum(0.5, 1.0, 0.5, 3.9, 6.1)   # different support, same total mass
# match centroid explicitly for illustration of 'matched low moments':
print(f"   A  M0,M1 = {moments(xa,wa,2)}   B M0,M1 = {moments(xb,wb,2)}")
print(f"   lambda_2(1):  A={christoffel(xa,wa,1.0,2):.6f}   B={christoffel(xb,wb,1.0,2):.6f}")
print(f"   (lambda_2 depends only on M_0,M_1,M_2 -> matching those forces equality;")
print(f"    a thin atom can be absorbed into matched moments = SLOT-BOOTSTRAP Delta_min=0.)")

print("\n(D) ACID: atomless continuum from s=0.3 to 6 (no atom) -> lambda_n(1) must -> 0?")
xg, wg = make_atom_continuum(0.0, 1.0, 1.0, 0.3, 6.0)
for n in range(1, 9):
    print(f"   n={n:2d}: lambda_n(1) = {christoffel(xg, wg, 1.0, n):.6e}")

print("\n(E) is the Gauss w_1 a usable LOWER bound forcing f1^2>0?")
# Gauss-2 rule of A (matches M_0..M_3). That SAME rule is the 2-node measure with matched
# moments; a pure-continuum measure can also match those moments -> w_1>0 does NOT force atom.
nodes, wts = gauss(xa, wa, 2)
print(f"   A Gauss-2: nodes={nodes}  weights={wts}  (w_1={min(wts):.5f} ~ atom, lower)")
nodes2, wts2 = gauss(xg, wg, 2)
print(f"   atomless Gauss-2: nodes={nodes2}  weights={wts2}")
print("   -> the atomless measure ALSO yields a positive smallest weight: w_1>0 does NOT force an atom.")
print("=" * 74)
