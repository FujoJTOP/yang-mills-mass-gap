#!/usr/bin/env python3
r"""
tr_maxplus.py -- TR (tropical / min-plus / Maslov) species: does the max-plus eigenvalue
(cycle mean) of the coarse transfer kernel see the floor?

OBJECT (exact, reused from P4/D3/D1): the periodic constant-mode transfer kernel
    M_b[i,j] = sqrt(r1w_i r1w_j) exp(beta cos((z_j - z_i)/2)),   r1w = rho1(z)*quadw,
rho1(z) = pi (Cin(2pi) - Cin|z|) the EXACT one-link SU(2) Cartan Haar marginal.
Anchors to reproduce (P4-PERIODIC-TRANSFER): beta=0  l1/l0 = 5.7e-16 (rank 1);
beta=2.5 l1/l0 = 0.52498, Delta = -log(l1/l0) = 0.6444; Haar floor 1/Var(xi3) = 0.26883.

TROPICAL DATA.  W = (1/beta) log M  (so M = e^{beta W}).  With rho_inf = max cycle mean of W:
    [T]  e^{beta rho_inf} <= rho(M) <= e^{beta rho_inf} * N_eff ,     N_eff = max_i sum_j e^{beta(W_ij + u_j - u_i - rho_inf)}
(the Collatz-Wielandt / log-sum-exp = Maslov-dequantization squeeze, u = tropical eigenvector).
So (1/beta) log lambda_0 in [rho_inf, rho_inf + (log N_eff)/beta]: a STRICT two-sided bound on the
DOMINANT eigenvalue only.  Direction account (charter family 14): m = log(l0/l1) needs an UPPER
bound on lambda_1 and a LOWER bound on lambda_0.  The squeeze supplies the lambda_0 half only; the
lambda_1 half is supplied (if at all) by the Dobrushin ergodicity coefficient delta(P) >= l1/l0,
whose tropical limit is 1 whenever the ground-state map x -> argmax_y W(x,y) is non-constant.

SECTIONS
  [A] anchors: beta=0 rank-1, beta=2.5 Delta=0.6444, Haar floor 0.26883.
  [B] max-plus: rho_inf (Karp), critical-graph degeneracy, exact N_eff, squeeze width vs Delta.
  [C] degeneracy of the critical graph at beta -> infinity: W -> cos((z-z')/2), all self-loops tie
      => rho_inf = 1 with an n-fold degenerate maximiser, tropical Delta = 0.
  [D] grid refinement n: does the squeeze width shrink (good) or grow (fatal)?
  [E] Dobrushin delta(P) vs the exact l1/l0: the only tropical-compatible UPPER bound on lambda_1.
  [F] beta-scan of the exact Delta: the tropical limit (beta -> infinity) is the CONTINUUM limit,
      where the lattice-units gap closes => the mechanism's own limit is the anti-target.

Run: py -3.12 -I tr_maxplus.py [quick|full].  numpy/scipy only, deterministic, no MC.
Labels [T]/[V].  No QED.  Clay 0%.
"""
import sys
import numpy as np
from scipy.special import sici

np.seterr(invalid='ignore', divide='ignore')
EULER = 0.5772156649015328606
TWO_PI = 2.0 * np.pi


def Cin(x):
    x = np.asarray(x, float)
    _, ci = sici(x)
    out = EULER + np.log(np.where(x > 0, x, 1.0)) - ci
    return np.where(x > 0, out, 0.0)


def rho1(z):
    z = np.asarray(z, float)
    az = np.abs(z)
    return np.where(az <= TWO_PI, np.pi * (Cin(TWO_PI) - Cin(np.minimum(az, TWO_PI))), 0.0)


def nodes(nz):
    xg, wg = np.polynomial.legendre.leggauss(nz)
    z = TWO_PI * xg
    wz = TWO_PI * wg
    return z, rho1(z) * wz


def kernel(z, r1w, beta):
    s = np.sqrt(np.maximum(r1w, 0.0))
    K = np.exp(beta * np.cos(0.5 * (z[None, :] - z[:, None])))
    M = (s[:, None] * K) * s[None, :]
    return 0.5 * (M + M.T)


# ------------------------------------------------------------------ max-plus tools
def karp_cycle_mean(W):
    """Karp's algorithm: rho_inf = max over cycles of the mean weight. O(n^3), vectorised."""
    n = W.shape[0]
    F = np.zeros(n)                      # F_k(i): best weight of a walk of length k ending at i
    best = np.full(n, np.inf)            # min_k (F_n(i) - F_k(i)) / (n-k)
    for k in range(n):
        best = np.minimum(best, (F - F) * 0)  # placeholder, filled below
        break
    F = np.zeros(n)
    Fs = np.empty((n + 1, n))
    Fs[0] = F
    for k in range(n):
        # F_{k+1}(j) = max_i (F_k(i) + W_ij)
        Fs[k + 1] = np.max(Fs[k][:, None] + W, axis=0)
    Fn = Fs[n]
    ks = np.arange(n)
    with np.errstate(invalid='ignore'):
        vals = (Fn[None, :] - Fs[:n, :]) / (n - ks)[:, None]
    return float(np.max(np.min(vals, axis=0))), Fs


def tropical_eigvec(W, rho, iters=4000):
    """u solving max_j (W_ij + u_j) = rho + u_i (value iteration, normalised)."""
    n = W.shape[0]
    u = np.zeros(n)
    for _ in range(iters):
        un = np.max(W + u[None, :], axis=1) - rho
        un -= un.max()
        if np.max(np.abs(un - u)) < 1e-13:
            u = un
            break
        u = un
    resid = np.max(np.max(W + u[None, :], axis=1) - rho - u)
    return u, resid


def squeeze_width(M, W, rho, u):
    """exact Collatz-Wielandt / log-sum-exp width per row: (1/beta) log sum_j e^{beta(W_ij+u_j-u_i-rho)}."""
    beta = 1.0 / 1.0  # placeholder, overwritten by caller
    raise RuntimeError


def cw_width(M, W, rho, u, beta):
    """N_eff(i) = sum_j e^{beta(W_ij + u_j - u_i - rho)}; width_i = log N_eff(i)/beta.
    Return (max width, min width, N_eff array)."""
    arg = beta * (W + u[None, :] - u[:, None] - rho)
    mx = arg.max(axis=1)
    Neff = np.exp(arg - mx[:, None]).sum(axis=1) * np.exp(mx)
    w = np.log(Neff) / beta
    return float(w.max()), float(w.min()), Neff


def dobrushin(M):
    """delta(P) = (1/2) max_{i,i'} sum_j |P_ij - P_i'j|, P = row-normalised M.  >= l1/l0 for all."""
    rs = M.sum(axis=1)
    P = M / rs[:, None]
    n = P.shape[0]
    best = 0.0
    for i in range(n):
        d = np.abs(P - P[i][None, :]).sum(axis=1)
        best = max(best, d.max())
    return 0.5 * best


# ------------------------------------------------------------------ main
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'quick'
    nz = 320 if mode != 'full' else 460
    print("#" * 100)
    print("# TR -- tropical / max-plus species: cycle means vs the exact transfer spectrum (P4 kernel)")
    print("#" * 100)

    # exact one-link Haar variance (P4 anchor)
    zg = np.linspace(-TWO_PI, TWO_PI, 200001)
    wg = np.full(zg.size, 2 * TWO_PI / (zg.size - 1)); wg[[0, -1]] *= 0.5
    r1g = rho1(zg) * wg
    var1 = (zg * zg * r1g).sum() / r1g.sum()
    print(f"[A] exact one-link Var(xi3)_Haar = {var1:.6f}  => Haar floor 1/Var = {1/var1:.5f}")

    z, r1w = nodes(nz)
    print(f"[A] grid nz = {nz};  r1w: max={r1w.max():.6e}  sum={r1w.sum():.6f}")

    # ---------------- [A] anchors
    print("\n" + "=" * 100)
    print("[A] ANCHORS on the P4 kernel (must reproduce 5.7e-16 / 0.6444)")
    print("-" * 100)
    for beta in (0.0, 2.5):
        M = kernel(z, r1w, beta)
        ev = np.linalg.eigvalsh(M)[::-1]
        ratio = ev[1] / ev[0] if ev[0] > 0 else 0.0
        D = -np.log(max(ratio, 1e-300))
        print(f"    beta={beta:4.2f}: lam0={ev[0]:.6f} lam1/lam0={ratio:.6e} Delta={D:.4f} "
              f"top4={np.array2string(ev[:4], precision=4)}")

    # ---------------- [B] tropical data at beta=2.5
    print("\n" + "=" * 100)
    print("[B] MAX-PLUS DATA at beta=2.5: rho_inf (Karp), degeneracy, exact squeeze width vs Delta")
    print("-" * 100)
    beta = 2.5
    M = kernel(z, r1w, beta)
    ev = np.linalg.eigvalsh(M)[::-1]
    D_ex = -np.log(ev[1] / ev[0])
    W = np.log(M) / beta
    rho, Fs = karp_cycle_mean(W)
    u, resid = tropical_eigvec(W, rho)
    wmax, wmin, Neff = cw_width(M, W, rho, u, beta)
    # degeneracy of the critical graph: number of 1-cycles (self loops) attaining rho
    diag = np.diag(W)
    ndeg = int(np.sum(np.abs(diag - rho) < 1e-9))
    print(f"    rho_inf (max cycle mean)      = {rho:.6f}   [Karp, resid of eigenvector {resid:.2e}]")
    print(f"    tropical Delta_pred          = 0 exactly (critical graph: {ndeg} of {nz} self-loops tie at rho)")
    print(f"    exact  (1/beta)log lam0      = {np.log(ev[0])/beta:.6f}")
    print(f"    squeeze: [(1/b)log lam0] in [rho_inf, rho_inf + width],  width = log(N_eff)/beta")
    print(f"       N_eff: max={Neff.max():.4f}  min={Neff.min():.4f}   width: max={wmax:.4f} min={wmin:.4f}")
    print(f"    EXACT Delta (signal)         = {D_ex:.4f}")
    print(f"    ==> width/signal              = {wmax/D_ex:.3f}   (1.0 = the bound is as wide as the whole gap)")
    print(f"    lambda_0 bracket: [{np.exp(beta*rho):.4f}, {np.exp(beta*(rho+wmax)):.4f}]  contains {ev[0]:.4f}"
          f" -> {np.exp(beta*rho) <= ev[0] <= np.exp(beta*(rho+wmax))}")

    # ---------------- [C] beta -> infinity: the tropical limit is degenerate
    print("\n" + "=" * 100)
    print("[C] beta -> oo: W -> cos((z-z')/2); S_W is EXACTLY FLAT on the constant mode (P4 sec 1)")
    print("    => every self-loop ties at rho_inf; the tropical spectrum is n-fold degenerate;")
    print("    tropical Delta = 0 for EVERY beta.  Exact Delta(beta) (lattice units) below:")
    print("-" * 100)
    print(f"    {'beta':>6}{'lam0':>12}{'lam1/lam0':>12}{'Delta_exact':>13}{'rho_inf':>10}{'width':>9}{'width/Delta':>13}")
    for bb in (0.0, 0.5, 1.0, 2.5, 5.0, 10.0, 20.0):
        Mb = kernel(z, r1w, bb)
        evb = np.linalg.eigvalsh(Mb)[::-1]
        Db = -np.log(max(evb[1] / evb[0], 1e-300))
        if bb > 0:
            Wb = np.log(Mb) / bb
            rb, _ = karp_cycle_mean(Wb)
            ub, _ = tropical_eigvec(Wb, rb)
            wb, _, _ = cw_width(Mb, Wb, rb, ub, bb)
        else:
            rb, wb = np.nan, np.nan
        print(f"    {bb:6.2f}{evb[0]:12.5f}{evb[1]/evb[0]:12.5e}{Db:13.4f}{rb:10.5f}{wb:9.5f}"
              f"{(wb/Db if Db > 1e-12 else np.inf):13.4f}")

    # ---------------- [D] grid refinement
    print("\n" + "=" * 100)
    print("[D] GRID REFINEMENT at beta=2.5: does the squeeze width shrink (uniform) or grow (fatal)?")
    print("-" * 100)
    print(f"    {'nz':>5}{'lam0':>12}{'Delta':>9}{'rho_inf':>10}{'N_eff max':>11}{'width':>9}"
          f"{'width/Delta':>13}{'log(nz)/beta':>14}")
    for nzz in (40, 80, 160, 320, 640):
        zz, rr = nodes(nzz)
        Mz = kernel(zz, rr, 2.5)
        evz = np.linalg.eigvalsh(Mz)[::-1]
        Dz = -np.log(evz[1] / evz[0])
        Wz = np.log(Mz) / 2.5
        rz, _ = karp_cycle_mean(Wz)
        uz, _ = tropical_eigvec(Wz, rz)
        wz, _, Nz = cw_width(Mz, Wz, rz, uz, 2.5)
        print(f"    {nzz:5d}{evz[0]:12.5f}{Dz:9.4f}{rz:10.5f}{Nz.max():11.3f}{wz:9.4f}"
              f"{wz/Dz:13.3f}{np.log(nzz)/2.5:14.3f}")

    # ---------------- [E] Dobrushin: the only tropical-compatible UPPER bound on lambda_1
    print("\n" + "=" * 100)
    print("[E] DOBRUSHIN delta(P) >= lam1/lam0  (the only upper bound on lambda_1 in the tropical frame)")
    print("-" * 100)
    print(f"    {'beta':>6}{'lam1/lam0':>12}{'delta(P)':>11}{'-log delta':>12}{'exact Delta':>13}{'bound gap':>11}")
    for bb in (0.5, 1.0, 2.5, 5.0):
        Mb = kernel(z, r1w, bb)
        evb = np.linalg.eigvalsh(Mb)[::-1]
        ratio = evb[1] / evb[0]
        Db = -np.log(ratio)
        dl = dobrushin(Mb)
        print(f"    {bb:6.2f}{ratio:12.5f}{dl:11.5f}{-np.log(dl):12.4f}{Db:13.4f}"
              f"{(-np.log(dl)/Db if Db > 1e-12 else np.inf):11.3f}")
    print("    (bound gap = (-log delta)/Delta >= 1: how many times weaker the rigorous bound is)")

    print("\nVERDICT: dump to TR-TROPICAL.md.  No QED.  Clay 0%.")


if __name__ == "__main__":
    main()
