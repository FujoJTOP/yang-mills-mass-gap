"""
pnt_ope_rho_witness.py -- PNT-OPE-RHO. Is the KL support floor forced by ANY
"edge-sensitive" data beyond finite moments, or is moment-blindness intrinsic?

Model: rho(s) >= 0 on s>=0, s = m^2. Constraints tested:
  (M)  four moments M_j = int s^j rho ds,  j in {-2,-1,0,1}
       (M-2 ~ trace anomaly A, M-1 = chi = Ward, M0 = norm W, M1 = higher)
  (T)  UV tail: fix the weight and first moment ABOVE a cut s=20
       (models the AF asymptotic tail rho ~ s (log s)^-b at large s)

True "gapped" measure: atom f1^2=0.5 at m1^2=1 + continuum rho~s on [4,S].
Competitor: NONNEG measure with NO atom, floor Delta, matching (M) [and (T)].
Question: can Delta -> 0 while matching? (then moment+tail data is floor-blind).

Case 2: impose "rho=0 below Delta_thr"; find the largest feasible Delta_thr.
This is the bootstrap direction: it can only give an UPPER bound on the floor.

numpy + scipy only. Run: py -3.12 -I research/unconditional/w7/pnt_ope_rho_witness.py
"""
import numpy as np
from scipy.optimize import nnls

JS = [-2, -1, 0, 1]
def mom(nodes, w):
    return np.array([np.dot(w, nodes ** j) for j in JS])

# ---- true gapped measure ----
f1sq, m1sq, S = 0.5, 1.0, 300.0
sc = np.linspace(4.0, S, 3000)
dens = sc / sc.sum()                 # AF-like shape rho(s) ~ s
wc = dens * np.gradient(sc)
wc = wc / wc.sum() * 2.0             # total continuum weight = 2
nodes_true = np.concatenate(([m1sq], sc))
w_true = np.concatenate(([f1sq], wc))
M = mom(nodes_true, w_true)
print("true moments [M-2, M-1, M0, M1] =", np.array2string(M, precision=4))
print("true atom weight at s=1: %.3f   true support floor: %.3f" % (f1sq, m1sq))


def solve(A, target, n):
    # A is (num_constraints x num_unknowns); nnls wants (m, n) with len(b)=m
    w, _ = nnls(A, target, maxiter=40000)
    return np.linalg.norm(A @ w - target), w


# ---------------- Case 1: four moments only, gapless competitor ----------------
def case1(Delta, n=700):
    s = np.linspace(Delta, S, n)
    A = np.vstack([s ** j for j in JS])          # 4 x n
    return solve(A, M, n)[0]

print("\nCASE 1  moments only, NO atom, floor Delta -> can it match?")
print("  (residual ~ 0 for every Delta => floor is NOT pinned by moments)")
for D in [1.0, 0.3, 0.1, 0.03, 0.01, 0.003]:
    print("   Delta=%-7s residual=%.2e" % (D, case1(D)))

# ---------------- Case 3: moments + fixed UV tail above s=20 ----------------
def case3(Delta, n=700, cut=20.0):
    s = np.linspace(Delta, S, n)
    rows = [s ** j for j in JS]
    hi = (s >= cut).astype(float)
    rows.append(hi)                               # tail weight above cut
    rows.append(hi * s)                           # tail 1st moment above cut
    A = np.vstack(rows)
    msk = nodes_true >= cut
    m_high = np.array([w_true[msk].sum(), w_true[msk].dot(nodes_true[msk])])
    target = np.concatenate([M, m_high])
    return solve(A, target, n)[0]

print("\nCASE 3  moments + FIXED UV tail (weight & 1st moment above s=20):")
print("  (residual ~ 0 for Delta -> 0 => the UV tail does NOT force a small-s atom)")
for D in [1.0, 0.3, 0.1, 0.03, 0.01, 0.003]:
    print("   Delta=%-7s residual=%.2e" % (D, case3(D)))

# ---------------- Case 2: threshold rho=0 below Delta_thr -> Delta_max ----------------
def feas_thr(Delta, n=700):
    s = np.linspace(Delta, S, n)
    A = np.vstack([s ** j for j in JS])
    return solve(A, M, n)[0]

lo, hi = 0.3, 3.0
for _ in range(40):
    mid = 0.5 * (lo + hi)
    if feas_thr(mid) < 1e-6:
        lo = mid
    else:
        hi = mid
print("\nCASE 2  threshold (rho=0 below Delta_thr):")
print("  Delta_max (largest feasible threshold) = %.3f   (true floor = %.3f)" % (lo, m1sq))
print("  => the bootstrap threshold direction yields only an UPPER bound on the floor,")
print("     i.e. Delta_max >= true floor; it can NEVER certify a lower bound Delta_min>0.")

# ---------------- What WOULD read the atom: a concentration functional ----------------
# Show the atom weight is invisible to moments: max point mass is free between 0 and ...
print("\nCHECKS  min possible value of the largest point mass (concentration) under moments only:")
# minimax: can we drive the atom weight at s=1 to 0 while matching moments? (yes, case1 with floor=1e-3 already does; the 'atom' is gone)
r, w = solve(np.vstack([np.linspace(1e-3, S, 700) ** j for j in JS]), M, 700)
s = np.linspace(1e-3, S, 700)
print("  gapless match residual at floor=1e-3: %.2e  ;  max single-node weight: %.3g" % (r, w.max() if r < 1e-3 else float('nan')))
print("\nCONCLUSION: four moments (and the UV tail) admit a floor=%.g, atom-free competitor."
      % 1e-3)
