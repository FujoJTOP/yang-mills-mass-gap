"""Ensemble law-of-total-variance proxy for the Var term of the marginal Hessian.

G''(A') = E[P^T H P | A'] - Var(P^T dS | A').

Exact constrained sampling given A' is infeasible on SU(3) without a constrained
MCMC.  This script computes the ENSEMBLE-average estimate:

  E_A'[G'' ] = E[P^T H P] - E_A'[ Var(P^T dS | A') ]
             = Hbar - ( Var(g) - Var(E[g|A']) )          (law of total variance)
  with  g = P^T dS  (coarse gradient),  A' = P^T xi  (coarse field).

Var(E[g|A']) is estimated by a LINEAR regression of g on A' (a lower bound on the
true conditional-mean variance -> an UPPER bound on the Var term).  We do the
regression in the SOFTest eigen-direction v of Hbar (1-D, well-posed with ~60
configs; the full Nc=512 covariance is rank-deficient).

Coordinate caveat: A' uses the ABSOLUTE principal log of the link; g uses the
left-variation gradient.  These are mixed coordinates, so this is a PROXY.
"""
import sys

import numpy as np

sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, 'D:/Dev/math/research/cpb-attack/11-slc-blocked')
import su3 as S  # noqa: E402
from full_hessian import build_hessian  # noqa: E402
from coarse_block_hessian import lift_matrix  # noqa: E402


def log_xi(U):
    """Absolute algebra coords: U = exp(i xi^a T^a), xi^a = Im tr(T^a log U)/2."""
    w, V = np.linalg.eig(U)
    th = np.angle(w)                      # (...,3)
    Vin = np.linalg.inv(V)
    logU = np.einsum('...ij,...j,...jk->...ik', V, 1j * th, Vin)
    return 0.5 * np.einsum('aij,...ji->...a', S.T_GEN, logU).imag


def main():
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    ncfg = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    nH = int(sys.argv[3]) if len(sys.argv) > 3 else 15
    beta = 2.5
    rng = np.random.default_rng(31)
    lat = S.Lattice(L, beta, rng, cold=True)
    for _ in range(300):
        lat.sweep()
    P = lift_matrix(L, 2)
    P = P.tocsr()

    As, Gs = [], []
    Hbar = None
    for c in range(ncfg):
        lat.sweep()
        U = lat.U
        Aflat = log_xi(U).reshape(-1)
        gflat = S.action_gradient(U, beta).reshape(-1)
        As.append(P.T @ Aflat)
        Gs.append(P.T @ gflat)
        if c < nH:
            H = build_hessian(U, beta, L)
            Hc = (P.T @ (H @ P)).toarray()
            Hc = 0.5 * (Hc + Hc.T)
            Hbar = Hc if Hbar is None else Hbar + Hc
        if (c + 1) % 10 == 0:
            print(f"  [{c+1}/{ncfg}]", flush=True)

    A = np.array(As)                      # (ncfg, Nc)
    G = np.array(Gs)                      # (ncfg, Nc)
    Hbar /= nH
    ev, evec = np.linalg.eigh(Hbar)
    v = evec[:, 0]                        # softest E-direction
    lam_E = ev[0]
    print(f"\nHbar (n={nH}) dim={ev.size}: min={lam_E:+.4f}  "
          f"lowest5={np.round(ev[:5],4)}")

    gv = G @ v
    Av = A @ v
    # 1-D linear regression gv = a + b Av, residual variance
    Avc = Av - Av.mean()
    gvc = gv - gv.mean()
    b = (Avc @ gvc) / (Avc @ Avc) if (Avc @ Avc) > 0 else 0.0
    resid = gvc - b * Avc
    var_g = float(np.var(gv))
    var_resid = float(np.var(resid))
    var_expl = var_g - var_resid
    R2 = var_expl / var_g if var_g > 0 else 0.0
    Gpp_soft = lam_E - var_resid

    print(f"  softest direction v:  lam_E = {lam_E:+.4f}")
    print(f"  Var(gv)={var_g:.4f}  residual Var={var_resid:.4f}  "
          f"explained={var_expl:.4f}  R^2={R2:.4f}")
    print(f"  G''_soft = lam_E - residualVar = {lam_E:+.4f} - {var_resid:.4f} "
          f"= {Gpp_soft:+.4f}")
    print(f"  (earlier within-block/16 proxy ~ 0.19)")
    verdict = "DEAD-lean (G''_soft<0)" if Gpp_soft < 0 else "ALIVE-lean (G''_soft>0)"
    print(f"\nVERDICT: {verdict}; but residual Var is an UPPER bound "
          f"(linear model), mixed coordinates -> proxy only.")

    # full-covariance diagnostics (rank-deficient descriptive only)
    try:
        Cg = np.cov(G.T)
        CA = np.cov(A.T)
        print(f"\n  descriptive: max eig Cov(g)={np.linalg.eigvalsh(Cg)[-1]:.3f}, "
              f"max eig Cov(A')={np.linalg.eigvalsh(CA)[-1]:.3f}")
    except Exception as e:  # noqa: BLE001
        print("  covariance diag skipped:", e)


if __name__ == '__main__':
    main()
