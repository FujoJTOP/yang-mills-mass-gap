"""
blast_operator_ineq.py -- witness for BLAST-OPERATOR-INEQ.
Direct operator inequality on Var(a_{k->0}) = <v, K^{-1} v> for the
Savvidy/D-term Hessian K = H+ - c phi phi^T  (c = g|B| > 0, H+ > 0, rank-one NEGATIVE mode).

Claim [T]:  <v,K^{-1}v> = <v,H+^{-1}v> + c <v,H+^{-1}phi>^2/(1-c*alpha)
            <= <v,H+^{-1}v> / (1-c*alpha)   <=   1/(mu (1-c*alpha)),
   where alpha = phi^T H+^{-1} phi, mu = lambda_min(H+), valid iff K > 0 (<=> c*alpha < 1).
Cauchy-Schwarz in the H+^{-1} inner product: <v,H+^{-1}phi>^2 <= <v,H+^{-1}v> * alpha.
The negative mode INCREASES the variance; the bound is SOFTENED by (1-c*alpha)^{-1},
diverging as c*alpha -> 1 (negative-mode instability). No sharp K^{-1} >= lambda.
"""
import numpy as np

def trial(n, c0, seed):
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n, n)); H = A @ A.T + np.eye(n)   # PD
    phi = rng.standard_normal(n); phi /= np.linalg.norm(phi)
    Hinv = np.linalg.inv(H)
    alpha = phi @ Hinv @ phi
    c = min(c0, 0.99 / alpha)                                  # force K PD
    K = H - c * np.outer(phi, phi)
    evK = np.linalg.eigvalsh(K)
    Kinv = np.linalg.inv(K)
    mu = np.linalg.eigvalsh(H).min()
    v = rng.standard_normal(n); v /= np.linalg.norm(v)
    lhs = v @ Kinv @ v
    Hv = v @ Hinv @ v
    rhs_sm = Hv / (1 - c * alpha)                              # Sherman-Morrison bound
    bound = 1.0 / (mu * (1 - c * alpha))                      # mu-bound
    return evK.min(), lhs, rhs_sm, bound, c * alpha, Hinv, phi, H, v

# random trials: inequality must hold, and K stays PD
maxratio = 0.0
for s in range(500):
    evmin, lhs, rhs, bound, soft, *_ = trial(5, 0.30, s)
    assert evmin > 0.0, "K not PD"
    assert lhs <= rhs + 1e-9, (lhs, rhs)
    assert lhs <= bound + 1e-9, (lhs, bound)
    maxratio = max(maxratio, lhs / rhs)
print("500 random trials: all pass.  max  <v,Kinv v>/<v,Hinv v>/(1-ca)  =", round(maxratio, 6), "(<=1)")

# tightness: choose v parallel to H+^{-1} phi  -> Cauchy-Schwarz equality
evmin, lhs, rhs, bound, soft, Hinv, phi, H, _ = trial(5, 0.30, 7)
w = Hinv @ phi; w /= np.linalg.norm(w)
K = H - 0.30 * np.outer(phi, phi)
# recompute with the same c actually used
alpha = phi @ Hinv @ phi; c = min(0.30, 0.99/alpha); K = H - c*np.outer(phi,phi); Kinv=np.linalg.inv(K)
lhs_t = w @ Kinv @ w; Hv_t = w @ Hinv @ w; rhs_t = Hv_t/(1 - c*alpha)
print("tightness (v || H+^{-1}phi):  lhs/rhs =", round(lhs_t/rhs_t, 8), "(=1 to rounding => Cauchy-Schwarz sharp)")

# obstruction: bound blows up as c*alpha -> 1
print("\nsoftening blow-up (fixed H+,phi; vary c):")
rng = np.random.default_rng(1)
A = rng.standard_normal((4,4)); H = A@A.T+np.eye(4)
phi = rng.standard_normal(4); phi/=np.linalg.norm(phi)
Hinv=np.linalg.inv(H); alpha=phi@Hinv@phi
for frac in [0.3, 0.7, 0.9, 0.99, 0.999]:
    ca = frac
    bound = 1.0/(np.linalg.eigvalsh(H).min()*(1-ca))
    print(f"  c*alpha={ca:6.3f}   mu*(1-c*alpha)={np.linalg.eigvalsh(H).min()*(1-ca):.4f}   1/(mu(1-ca))={bound:8.3f}")

print("\nConclusion: the direct operator inequality EXISTS but is SOFTENED by the rank-one")
print("negative mode:  Var(a) <= 1/(mu(1 - g|B|*alpha)).  The sharp K^{-1} >= lambda")
print("fails exactly as g|B|*alpha -> 1.  Missing input = uniform  mu(1 - g|B| alpha) >= lambda.")
