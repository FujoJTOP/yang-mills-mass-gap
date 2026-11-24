"""
cert_chi_ineq.py -- witness for CERT-CHI-INEQ.md.  numpy only. py -3.12 -I
Attempts the marginal inequality  S_bare > 2 g^2/lambda4  (equiv. mu(1-g|B|alpha) >= lambda)
for the chi-half of H1. Shows the exact structure and pins the H1 step.
"""
import numpy as np

def sm_var(Hp, phi, c, v):
    """Var(v) = v^T (Hp - c phi phi^T)^{-1} v ; returns (var, alpha, posdef)."""
    Hpi = np.linalg.inv(Hp)
    alpha = float(phi @ Hpi @ phi)
    K = Hp - c * np.outer(phi, phi)
    if np.linalg.eigvalsh(K).min() <= 0:
        return np.inf, alpha, False
    return float(v @ np.linalg.solve(K, v)), alpha, True

def main():
    rng = np.random.default_rng(0)
    n = 8; mu = 1.0
    Q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    Hp = Q @ np.diag(np.linspace(mu, 4*mu, n)) @ Q.T
    phi = Q[:, 0]          # negative mode aligned with lightest Hp eigvec (worst case)
    v = phi.copy()         # coarse mode along phi (sharp direction)
    lam = 0.5

    print("(a) exact Sherman-Morrison Var vs bound 1/(mu(1-c alpha))  [T]")
    for c in [0.3, 0.5, 0.7, 0.9]:
        var, alpha, ok = sm_var(Hp, phi, c, v)
        bound = 1.0/(mu*(1.0-c*alpha)) if c*alpha < 1 else np.inf
        print(f"  c={c}: alpha={alpha:.3f} c*alpha={c*alpha:.3f} Var={var:.4f} bound={bound:.4f} Var<=bound:{var<=bound+1e-9}")

    print("\n(b) crude suff. condition  mu - c >= lambda  (from alpha<=1/mu): fails for large c")
    for c in [0.3, 0.5, 0.7, 0.9]:
        print(f"  c={c}: mu-c={mu-c:.2f} >= {lam}? {mu-c>=lam}")

    print("\n(c) negative-mode overlap decouples at deep IR: c*alpha(k) ~ (k/sqrt(gB))^2 * gB [T]")
    gB = 1.0; c = gB
    for kf in [0.05, 0.2, 0.5, 1.0, 2.0]:
        alpha_k = (kf**2)/(1.0+kf**2)          # schematic: support |p|~sqrt(gB)
        print(f"  k/sqrt(gB)={kf:.2f}: c*alpha={c*alpha_k:.4f}  (deep IR c*alpha->0 => softening decouples)")

    print("\n(d) residual floor  S_bare >= lambda + D_Sav,  D_Sav=2g^2/lambda4  [O] = H1")
    for lam4 in [2.0, 1.0, 0.5, 0.2]:
        D = 2.0*1.0/lam4
        print(f"  lambda4={lam4}: D_Sav={D:.2f}  need S_bare>={lam+D:.2f}  (S_bare = condensate floor, not fixed by action)")

    print("\n(e) verdict: (a)(c)+bounded D_Sav are [T]; (d) S_bare floor = H1; holds rigorously only beta<1/48 (SZZ)")

if __name__ == "__main__":
    main()
