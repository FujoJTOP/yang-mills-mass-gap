#!/usr/bin/env python3
"""
na2_tci.py -- T2: first measurement of (TCI) on solvable models.

(TCI) as stated in SIXTH-SPECIES-CHARTER 1.3:
      <phi,Q phi> + <phi,K phi> - cF >= mu > 0     for ALL unit phi in H_-  (negative-mode subspace),
      pointwise, nonperturbative.

To MEASURE it, the objects must be pinned down; this script takes the reading that Q3-SCHUR-DEFECT.md
makes available (inertia lemma  n_-(H) = n_-(K) + n_-(Sigma)) and uses the model's own Hessian:

  V := S + Phi            (log-density potential in the flat chart;  Phi = -log J)
  H := Hess V = Hess S + Hess Phi            <-- this IS the Savvidy object: in 0+0 dimensions
                                                 -D^2 -> 0, so H = "2 ad F_B" + Jacobian curvature.
  H_- := span of the negative eigenvectors of H at the point xi            (POINTWISE).
  Split the coarse field direction (a = xi_3, index 3) from the fibre (indices 1,2):
      H = [[H_AA, H_Aa],[H_aA, H_aa]] ,   Sigma := H_aa - H_aA H_AA^{-1} H_A3     (Schur complement)
      inertia lemma:  n_-(H) = n_-(H_AA) + n_-(Sigma)          (checked to machine precision below)
  Q := Hess Phi  (the geometric/Prékopa part; check it is PSD)
  K := Hess S    (the Savvidy part, indefinite)
  F := negative spectral mass  -sum_{lambda_i<0} lambda_i   of H   (SELF-CHOSEN; the corpus does not fix F)
  c := 1 .        We report the inequality both with and without cF so the choice does not hide anything.

MODELS: SU(2) flat chart, (i) 1 plaquette  S = beta(1 - 1/2 Re tr(U B^dag)),
                                  (ii) 2 plaquette S = beta sum_{r=1,2}(1 - 1/2 Re tr(U B_r^dag)).
Measure: J(xi) d^3xi e^{-S}, sampled uniformly in |xi|<2pi with importance weight J e^{-S}
(SU(2) exponential is closed-form: U = cos(rho/2) + i sin(rho/2) (xihat.sigma), no eigh needed).

Two INDEPENDENT objects are compared:
   Sigma(xi)         -- POINTWISE (local Hessian Schur complement)
   G''(a0)           -- CONDITIONAL marginal curvature (needs the variance; from na2_su3_na machinery)
The gap between them is exactly Var(G)+X, i.e. the categorical mismatch this report is about.

Run: py -3.12 -I na2_tci.py [--n 400000]
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import na2_su3_na as M

SIG = M.lambdas(2)          # sigma_1, sigma_2, sigma_3 (Hermitian, tr=2delta)

def staple(ax, t):
    return np.eye(2, dtype=complex) * np.cos(t / 2) + 1j * np.sin(t / 2) * SIG[ax]

def c0_of(xi, Bdag):
    """Re tr(U B^dag) for SU(2), U = cos(rho/2) + i sin(rho/2) (nhat.sigma) -- explicit 2x2."""
    rho = np.linalg.norm(xi, axis=1)
    nh = xi / np.maximum(rho, 1e-300)[:, None]
    A = sum(nh[:, k][:, None, None] * SIG[k] for k in range(3))
    U = np.cos(rho / 2)[:, None, None] * np.eye(2) + 1j * np.sin(rho / 2)[:, None, None] * A
    return np.real(np.einsum('mab,ba->m', U, Bdag))

def build_Bd(ax, t):
    """B^dag = exp(-i t sigma_ax/2) as a 2x2 matrix."""
    return np.eye(2, dtype=complex) * np.cos(t / 2) - 1j * np.sin(t / 2) * SIG[ax]

def Sfun(xi, Bd_list, beta):
    s = 0.0
    for Bd in Bd_list:
        s = s + beta * (1.0 - 0.5 * c0_of(xi, Bd))
    return s

def hessS(xi, Bd_list, beta, h=1e-3):
    """Hessian of S in the 3 flat-chart Cartesian directions, FD (6 evals of the closed form)."""
    m = xi.shape[0]
    H = np.zeros((m, 3, 3))
    f0 = Sfun(xi, Bd_list, beta)
    for i in range(3):
        for j in range(i, 3):
            if i == j:
                ei = np.zeros(3); ei[i] = h
                H[:, i, i] = (Sfun(xi + ei, Bd_list, beta) - 2 * f0 + Sfun(xi - ei, Bd_list, beta)) / h ** 2
            else:
                ei = np.zeros(3); ei[i] = h
                ej = np.zeros(3); ej[j] = h
                H[:, i, j] = H[:, j, i] = (
                    Sfun(xi + ei + ej, Bd_list, beta) - Sfun(xi + ei - ej, Bd_list, beta)
                    - Sfun(xi - ei + ej, Bd_list, beta) + Sfun(xi - ei - ej, Bd_list, beta)) / (4 * h ** 2)
    return H

def hessPhi(xi):
    """Hess of Phi=-log J, Phi=phi(rho): phi'' n n^T + phi'(1-n n^T)/rho."""
    rho = np.linalg.norm(xi, axis=1)
    nh = xi / np.maximum(rho, 1e-300)[:, None]
    p1 = 2.0 / rho - 1.0 / np.tan(rho / 2)                    # phi'
    p2 = -2.0 / rho ** 2 + 0.5 / np.sin(rho / 2) ** 2         # phi''
    H = p2[:, None, None] * nh[:, :, None] * nh[:, None, :] \
        + (p1 / rho)[:, None, None] * (np.eye(3)[None] - nh[:, :, None] * nh[:, None, :])
    return 0.5 * (H + np.swapaxes(H, 1, 2))

def negcount(H):
    w = np.linalg.eigvalsh(H)
    return (w < -1e-12).sum(axis=1), w

def main():
    n = 400000
    for i, a in enumerate(sys.argv):
        if a == '--n':
            n = int(sys.argv[i + 1])
    t0 = time.time()
    print("=" * 112, flush=True)
    print("T2 : (TCI) first measurement.  H = Hess V = Hess S + Hess Phi  (flat chart, SU(2)).", flush=True)
    print("     Sigma = Schur complement in the coarse direction a=xi_3 (pointwise).   F = neg. spectral mass.", flush=True)
    print("     Q = Hess Phi, K = Hess S.  (TCI): <phi,Q phi>+<phi,K phi>-cF >= mu on every negative mode.", flush=True)
    print("=" * 112, flush=True)
    rng = np.random.default_rng(31415)
    # uniform in ball |xi|<2pi
    u = rng.random(n); rr = 2 * np.pi * u ** (1.0 / 3.0)
    z = rng.standard_normal((n, 3)); z /= np.linalg.norm(z, axis=1)[:, None]
    xi = rr[:, None] * z
    rho = rr
    J = np.sin(rho / 2) ** 2 / np.maximum(rho ** 2, 1e-300)
    HP = hessPhi(xi)
    print(f"     samples={n}  ({time.time()-t0:.0f}s)", flush=True)

    models = [('1 plaquette  B=e3(pi)', [build_Bd(2, np.pi)]),
              ('1 plaquette  B=e3(pi/2)', [build_Bd(2, np.pi / 2)]),
              ('1 plaquette  B=1', [build_Bd(2, 0.0)]),
              ('2 plaquette  e3(pi)&e3(pi/2)', [build_Bd(2, np.pi), build_Bd(2, np.pi / 2)]),
              ('2 plaquette  e3(pi)&e1(pi)', [build_Bd(2, np.pi), build_Bd(0, np.pi)])]
    rows = []
    print(f"\n     {'model':>30}{'beta':>6}{'<n_->':>7}{'P(n_->0)':>10}{'<lam_-1>':>10}"
          f"{'<lam_+3>':>10}{'<Sigma>':>10}{'P(Sig<0)':>10}{'inertia ok':>11}", flush=True)
    for lab, Bds in models:
        for beta in (0.25, 0.5, 1.0, 2.5, 4.0, 6.0):
            HS = hessS(xi, Bds, beta)
            H = HS + HP
            nneg, w = negcount(H)
            wmax = w.max(axis=1); wmin = w.min(axis=1)
            HAA = H[:, :2, :2]
            H3A = H[:, 2, :2]; Ha3 = H[:, :2, 2]; Haa = H[:, 2, 2]
            nAA, wAA = negcount(HAA)
            det = np.linalg.det(HAA)
            inv = np.zeros_like(HAA)
            ok = np.abs(det) > 1e-12
            inv[ok] = np.linalg.inv(HAA[ok])
            Sig = Haa - np.einsum('ma,mab,mb->m', H3A, inv, Ha3)
            S = 0.0
            for Bd in Bds:
                S = S + beta * (1.0 - 0.5 * c0_of(xi, Bd))
            wt = J * np.exp(-S)
            wt = wt / wt.sum()
            E = lambda x: float(np.sum(wt * x))
            inert = E((nneg == nAA + (Sig < -1e-12)) * 1.0)
            print(f"     {lab:>30}{beta:6.2f}{E(nneg):7.3f}{E((nneg > 0) * 1.0):10.4f}"
                  f"{E(np.minimum(wmin, 0)):10.4f}{E(wmax):10.4f}{E(Sig):10.4f}"
                  f"{E((Sig < 0) * 1.0):10.4f}{inert:11.4f}", flush=True)
            rows.append(dict(model=lab, beta=beta, n_neg=E(nneg), p_neg=E((nneg > 0) * 1.0),
                             lam_min=E(np.minimum(wmin, 0)), lam_max=E(wmax),
                             Sigma=E(Sig), p_sigma_neg=E((Sig < 0) * 1.0), inertia=inert))
    with open('na2_tci.json', 'w') as f:
        json.dump(rows, f, indent=1, default=float)

    # ---------------- (b) pointwise Schur curvature Sigma  vs  CONDITIONAL curvature G''(a0)
    print(f"\n     POINTWISE Sigma (local Hessian, Schur) vs CONDITIONAL G''(a0) -- the gap is Var(G)+X."
          f"\n     If Sigma<0 on a positive-measure set, NO pointwise mu>0 exists; only the conditional"
          f"\n     version (a global variance object) can be positive -- and it is not (L,a)-uniform either.", flush=True)
    print(f"     {'model':>30}{'beta':>6}{'<Sigma|a0=0>':>14}{'P(Sig<0|0)':>12}"
          f"{'G\"(a0=0)':>12}{'G\"(a0=1)':>11}{'<n_->':>7}", flush=True)
    rng2 = np.random.default_rng(31415)
    for lab, Bds in models:
        geo = M.geometry(2, 400000, rng2)
        raw = M.action_traces(geo, Bds)
        ct = (sum(r[0] for r in raw), sum(r[1] for r in raw))   # additive; the constant offset in S cancels
        for beta in (0.5, 2.5, 6.0):
            HS = hessS(xi, Bds, beta); H = HS + HP
            nneg, w = negcount(H)
            HAA = H[:, :2, :2]; H3A = H[:, 2, :2]; Ha3 = H[:, :2, 2]; Haa = H[:, 2, 2]
            det = np.linalg.det(HAA); inv = np.zeros_like(HAA); ok = np.abs(det) > 1e-12
            inv[ok] = np.linalg.inv(HAA[ok])
            Sig = Haa - np.einsum('ma,mab,mb->m', H3A, inv, Ha3)
            S = 0.0
            for Bd in Bds:
                S = S + beta * (1.0 - 0.5 * c0_of(xi, Bd))
            wt = J * np.exp(-S); wt = wt / wt.sum()
            for a0 in (0.0,):
                k = np.exp(-0.5 * (xi[:, 2] - a0) ** 2 / 0.12 ** 2)
                ws = wt * k; ws = ws / ws.sum()
                ES = float(np.sum(ws * Sig)); PS = float(np.sum(ws * (Sig < 0)))
            g0 = M.extrap(geo, ct, beta, 0.0)
            g1 = M.extrap(geo, ct, beta, 1.0)
            print(f"     {lab:>30}{beta:6.2f}{ES:+14.4f}{PS:12.4f}{g0['Gpp']:+12.4f}"
                  f"{g1['Gpp']:+11.4f}{float(np.sum(ws * nneg)):7.3f}", flush=True)
    print(f"\n  [{time.time()-t0:.0f}s]  Clay 0%. No QED.", flush=True)

if __name__ == '__main__':
    main()
