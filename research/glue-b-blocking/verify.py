"""Independent verification for glue-b-blocking (GL2a criterion, GL2b direction).

Pure stdlib: Fraction for the exact criterion, a small Jacobi eigensolver for
the interlacing demonstrations.  Nothing here is a proof -- it is a check.
"""
from fractions import Fraction as F
import math

# ---------------------------------------------------------------- GL2a criterion
# 3-4-5 rotation R_ij, cos=3/5, sin=4/5 (sign convention FIXED here).
def R45(i, j, n=4, sign=+1):
    M = [[F(1) if a == b else F(0) for b in range(n)] for a in range(n)]
    c, s = F(3, 5), F(sign * 4, 5)
    M[i][i] = c; M[j][j] = c
    M[i][j] = -s; M[j][i] = s
    return M

def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]

def matT(A):
    return [list(r) for r in zip(*A)]

def diag(vals):
    n = len(vals)
    return [[F(vals[i]) if i == j else F(0) for j in range(n)] for i in range(n)]

def is_integral(A):
    return all(x.denominator == 1 for row in A for x in row)

def det4(A):
    # Laplace on 4x4, exact
    if len(A) == 1:
        return A[0][0]
    if len(A) == 2:
        return A[0][0]*A[1][1] - A[0][1]*A[1][0]
    tot = F(0)
    for j in range(len(A)):
        minor = [row[:j] + row[j+1:] for row in A[1:]]
        tot += (-1)**j * A[0][j] * det4(minor)
    return tot

def inverse(M):
    # 4x4 exact inverse via adjugate
    n = len(M)
    def minor_det(i, j):
        sub = [row[:j] + row[j+1:] for k, row in enumerate(M) if k != i]
        return det4(sub) if n == 4 else None
    cof = [[(-1)**(i+j) * minor_det(i, j) for j in range(n)] for i in range(n)]
    adj = matT(cof)
    d = det4(M)
    assert d != 0
    return [[adj[i][j] / d for j in range(n)] for i in range(n)]

def timed_reflection():
    return diag([-1, 1, 1, 1])   # x_0 -> -x_0

def spatial_reflection(axis):
    v = [1, 1, 1, 1]; v[axis] = -1
    return diag(v)

print("=" * 72)
print("GL2a: (B1) criterion  A in W cap R W R^-1  <=>  B^-1 A B integral")
print("=" * 72)
print("convention: R_ij has cos=3/5, sin=+4/5;  B = 5R ;  Lambda' = B Z^4")
print(f"{'block':8} {'plane':>10} {'B^TB=25I':>9} {'detB':>6} {'B integral':>11} "
      f"{'B^-1 A B integral':>18}  verdict")
for (i, j, name) in [(0, 1, "B_12"), (0, 2, "B_13"), (1, 3, "B_24")]:
    R = R45(i, j)
    B = [[5 * x for x in row] for row in R]
    ok_gram = all(matmul(matT(B), B)[a][b] == (25 if a == b else 0)
                  for a in range(4) for b in range(4))
    d = det4(B)
    for A, Alabel in [(timed_reflection(), "time"), (spatial_reflection(i), f"spatial e{i}")]:
        M = matmul(matmul(inverse(B), A), B)
        integral = is_integral(M)
        violating = [(a, b, M[a][b]) for a in range(4) for b in range(4)
                     if M[a][b].denominator != 1]
        verdict = "PASS" if integral else "FAIL"
        print(f"{name:8} {Alabel:>10} {str(ok_gram):>9} {str(d):>6} {str(is_integral(B)):>11} "
              f"{str(integral):>18}  {verdict}")
        if violating:
            print(f"         first non-integral entry: M[{violating[0][0]}][{violating[0][1]}]"
                  f" = {violating[0][2]}   (lcm of denominators = "
                  f"{math.lcm(*[x[2].denominator for x in violating])})")

# explicit vector check, as in the contract: A e_0 expressed in the block basis
print("\nexplicit: time reflection on first block vector b_0 = B e_0")
R = R45(0, 1); B = [[5 * x for x in row] for row in R]
b0 = [row[0] for row in B]
A = timed_reflection()
Ab0 = [sum(A[a][b] * b0[b] for b in range(4)) for a in range(4)]
# solve B alpha = A b0  (exact, 4x4 by hand: only first two coords nonzero)
import itertools
# brute force rational solve on the 2x2 corner
p, q = Ab0[0], Ab0[1]
b0x, b0y = b0[0], b0[1]
b1x, b1y = B[0][1], B[1][1]
d = b0x * b1y - b1x * b0y
alpha = (p * b1y - b1x * q) / d
beta = (b0x * q - p * b0y) / d
print(f"  b_0 = {b0}")
print(f"  A b_0 = {Ab0}")
print(f"  A b_0 = ({alpha}) b_0 + ({beta}) b_1      -> in Lambda'? {alpha.denominator==1 and beta.denominator==1}")

# ---------------------------------------------------------------- GL2b direction
def jacobi_eig(A, iters=100):
    """Symmetric eigen-decomposition by cyclic Jacobi (floats).  Returns sorted desc."""
    n = len(A)
    a = [[float(A[i][j]) for j in range(n)] for i in range(n)]
    for _ in range(iters):
        off = sum(a[i][j]**2 for i in range(n) for j in range(n) if i != j)
        if off < 1e-24:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-18:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1 if theta >= 0 else -1) / (abs(theta) + math.sqrt(theta*theta + 1))
                c = 1 / math.sqrt(t*t + 1); s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
    return sorted((a[i][i] for i in range(n)), reverse=True)

def compress(T, keep):
    """P T P for P = projection onto coordinate subspace 'keep' (list of indices)."""
    n = len(T)
    P = [[1.0 if (i == j and i in keep) else 0.0 for j in range(n)] for i in range(n)]
    PT = [[sum(P[i][k] * T[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return [[sum(PT[i][k] * P[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

print("\n" + "=" * 72)
print("GL2b: lambda_2(P T P) <= lambda_2(T)  (compression lowers lambda_2)")
print("=" * 72)
cases = [
    ("slow mode INSIDE block algebra", [[1,0,0,0],[0,.5,0,0],[0,0,.25,0],[0,0,0,.1]], [0,1,2], "equality"),
    ("slow mode OUTSIDE block algebra", [[1,0,0,0],[0,.5,0,0],[0,0,.25,0],[0,0,0,.1]], [0,2,3], "strict"),
    ("slowest vector mixed out",        [[1,0,0,0],[0,.5,0,0],[0,0,.25,0],[0,0,0,.1]], [0,3],   "strict"),
]
for label, T, keep, expect in cases:
    ev = jacobi_eig(T)
    cT = compress(T, keep)
    evc = [x for x in jacobi_eig(cT) if abs(x) > 1e-12]
    lam2, mu2 = ev[1], sorted(evc, reverse=True)[1] if len(evc) > 1 else 0.0
    print(f"{label:34} lambda_2={lam2:.6f}  mu_2={mu2:.6f}  "
          f"mu_2<=lambda_2? {mu2 <= lam2 + 1e-12}   [{expect}]")

# b = 2 time-blocking: Delta_c = -(1/(b a)) log mu_2^b = -(1/a) log mu_2  >= Delta_f
print("\nDelta_c >= Delta_f whenever mu_2 <= lambda_2  (b cancels):")
for (label, lam2, mu2) in [("inside", 0.5, 0.5), ("outside", 0.5, 0.25), ("mixed", 0.5, 0.1)]:
    Df = -math.log(lam2); Dc = -math.log(mu2)
    print(f"  {label:8} Delta_f={Df:.4f}  Delta_c={Dc:.4f}  Delta_c>=Delta_f? {Dc >= Df - 1e-12}")

# perturbation bound: |log mu - log mu'| <= eps/(mu-eps)
print("\ntruncation error bound |log mu - log mu'| <= eps/(mu-eps):")
for (mu, eps) in [(0.5, 0.01), (0.5, 0.1), (0.2, 0.05)]:
    mup = mu - eps
    print(f"  mu={mu}  eps={eps}  bound={eps/(mu-eps):.5f}  actual={abs(math.log(mu)-math.log(mup)):.5f}")

# ---------------------------------------------------------------- GL5 tensor check
print("\n" + "=" * 72)
print("GL5 kinematic: dim-4 H_hyp-scalars of F  (should be tr F^2 only)")
print("=" * 72)
pairs = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
# quadratic forms in F_{mu nu} invariant under all signed coordinate perms of 4 axes:
# a symmetric matrix Q on the 6-dimensional space of 2-forms; invariant iff
# Q[ab][cd] constant on orbits of the pairing -> only two orbits under S4 acting on pairs
# (disjoint pairs = 3, intersecting pairs = 12).  Check: any invariant Q is alpha*I + beta*(pairing)
# and the "disjoint vs intersecting" split is NOT S4-invariant as a matrix on pairs, so only I.
# Demonstrate: the candidate O(4)-breaking density  sum_i (F_{0i})^2 - (1/3) sum_{i<j}(F_{ij})^2
# is permuted to a different operator by the transposition (0 1) in S4.
def perm_pair(sigma, p):
    return tuple(sorted((sigma[p[0]], sigma[p[1]])))
sigma = {0:1, 1:0, 2:2, 3:3}   # swap axes 0 and 1
time_terms = {perm_pair(sigma, (0,i)) for i in (1,2,3)}
space_terms = {perm_pair(sigma, (i,j)) for i,j in pairs if i!=j and 0 not in (i,j)}
print(f"  under (0 1): {{(0,1),(0,2),(0,3)}} -> {sorted(time_terms)}")
print(f"               {{(1,2),(1,3),(2,3)}} -> {sorted(space_terms)}")
print("  -> the two families are exchanged: no S4-invariant combination with different"
      "\n     weights exists.  Only sum over ALL six pairs survives = tr F^2 (dim 4).")

# same statement with the REDUCED group that distinguishes time (S3 on space only):
sigma2 = {0:0, 1:1, 2:3, 3:2}
t2 = {perm_pair(sigma2, (0,i)) for i in (1,2,3)}
s2 = {perm_pair(sigma2, (i,j)) for i,j in pairs if i!=j and 0 not in (i,j)}
print(f"  with only spatial S3, swap (2 3): {sorted(t2)} / {sorted(s2)}  (fixed!)")
print("  -> anisotropic dim-4 operator IS allowed by (Z/2)^4 x S3: this is the"
      "\n     sharp falsification test (anisotropic lattice YM has a preferred frame).")
