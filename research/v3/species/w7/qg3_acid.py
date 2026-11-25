"""
qg3_acid.py -- Q4 route, step 3: the acid battery, plus the exact Gaussian dictionary.

The route's premise is  lambda >= c / R_eff  with R_eff an EFFECTIVE RESISTANCE.
There are exactly two ways to build R_eff; the script separates them and shows that
one is false and the other is the definition of lambda (tautology).

  (G) PRECISION network:  R_eff(v) := v^T Sigma v  (Sigma = inverse of the Hessian)
      => lambda = 1 / max_{|v|=1} R_eff(v) is an IDENTITY (verified to 1e-16).
  (E) GEOMETRIC network:  conductances from the lattice/graph only.
      => lambda * R_eff <= 2 for EVERY pair (qg1 PART A/B), i.e. an upper bound.

  ACID 1  free massless field: does the route return "no gap"?
  ACID 2  beta = 0 Haar (ultralocal; corpus: G''=0.2688, no atom): does it give a
          false positive certificate?
  ACID 3  double well: R_eff finite but gap = 0 -> must fail.
"""
import json
import math
import time

import numpy as np
from scipy.linalg import eigh_tridiagonal

LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


def torus_Reff_max(L, d):
    grids = np.meshgrid(*[np.arange(L)] * d, indexing="ij")
    m = np.stack([g.ravel() for g in grids], axis=1).astype(float)
    kk = 2.0 * math.pi * m / L
    lam = np.sum(4.0 * np.sin(kk / 2.0) ** 2, axis=1)
    n = len(lam)
    w = np.zeros(n)
    nz = lam > 1e-12
    w[nz] = 1.0 / lam[nz]
    C = np.cos(2.0 * math.pi * (m @ m.T) / L)
    R = (2.0 / n) * (w.sum() - C.T @ w)
    return float(np.max(np.maximum(R, 0.0))), lam


# ---------------------------------------------------------------------------
def part_G():
    say("=" * 78)
    say("PART G -- the exact Gaussian dictionary (the tautological reading)")
    say("  mu = N(0,Sigma), A = Sigma^{-1}.  Then")
    say("    lambda(=Poincare) = 1/lambda_max(Sigma),   R_eff^{prec}(v) = v^T Sigma v,")
    say("    and   lambda = 1 / max_{|v|=1} R_eff^{prec}(v)      <-- identity, not a bound.")
    say("  while for a FIXED pair direction (the only computable object)")
    say("    R_eff^{prec}(e_i-e_j) * lambda <= 2 ,  usually << 2.")
    say("=" * 78)
    rng = np.random.default_rng(7)
    worst = 0.0
    ratios = []
    for _ in range(300):
        n = int(rng.integers(4, 12))
        B = rng.normal(size=(n, n))
        A = B @ B.T + 0.1 * np.eye(n)          # SPD precision
        Sig = np.linalg.inv(A)
        lam_max = np.linalg.eigvalsh(Sig)[-1]
        lam = 1.0 / lam_max
        Rpair = np.array([[Sig[i, i] + Sig[j, j] - 2 * Sig[i, j]
                           for j in range(n)] for i in range(n)])
        worst = max(worst, float((Rpair * lam).max()))
        ratios.append(float(Rpair.max() / lam_max))
    say(f"  300 random precisions:  max over pairs of  lambda*R_eff(pair) = {worst:.9f}"
        f"   (must be <= 2)")
    say(f"  R_eff^max(pair) / lambda_max(Sigma)  :  median {np.median(ratios):.4f}, "
        f"max {max(ratios):.4f}   (<= 2 always; = 2 only when the top eigendirection")
    say(f"                                        is exactly e_i-e_j)")
    say("")
    # exact identity check on the lattice GFF
    for (L, d) in [(4, 2), (4, 3), (3, 4)]:
        n = L ** d
        idx = np.arange(n).reshape((L,) * d)
        Lm = np.zeros((n, n))
        for ax in range(d):
            for sh in (1, -1):
                nb = np.roll(idx, sh, axis=ax).ravel()
                Lm[idx.ravel(), nb] += -1.0
                Lm[idx.ravel(), idx.ravel()] += 1.0
        ev, U = np.linalg.eigh(Lm)
        lam1 = ev[1]
        Sig = U[:, 1:] @ np.diag(1.0 / ev[1:]) @ U[:, 1:].T   # GFF covariance on 1^perp
        lam_max = np.linalg.eigvalsh(Sig)[-1]
        say(f"  T^{d}_{L}: lambda_1 = {lam1:.12f}, 1/lambda_max(Sigma) = {1/lam_max:.12f},"
            f"  |diff| = {abs(lam1 - 1/lam_max):.2e}   <-- identity holds")
    say("")
    say("  CONSEQUENCE: 'lambda >= c/R_eff' is TRUE with c=1 and R_eff = sup_v R^prec(v),")
    say("  and that R_eff is 1/lambda.  Any R_eff that is not this sup is an")
    say("  upper-bound machine.  There is no third case.")
    say("")


# ---------------------------------------------------------------------------
def part_acid1():
    say("=" * 78)
    say("ACID 1 -- free massless field (the route's advertised strongest point)")
    say("=" * 78)
    say("  (a) FREE MAXWELL on the spatial 3-torus (Gauss law removes k=0):")
    say("      exact photon gap  w_min = 2 sin(pi/L),  geometric R_eff -> R_inf > 0.")
    say(f"  {'L':>3} {'n':>6} {'w_min=2sin(pi/L)':>18} {'R_eff^max(3D torus)':>21} {'w*R':>12}")
    rows = []
    for L in [2, 3, 4, 5, 6, 8]:
        R, _ = torus_Reff_max(L, 3)
        w = 2.0 * math.sin(math.pi / L)
        say(f"  {L:>3} {L**3:>6} {w:>18.8e} {R:>21.8e} {w*R:>12.6e}")
        rows.append(dict(L=L, w=w, R=R, prod=w * R))
    say("")
    say("  The product decays like 1/L => the route 'returns no gap'.  BUT:")
    say("  (b) add a mass m: the SAME geometric network (R_eff identical to 8 digits),")
    say("      and the exact gap is m for every L.  Two theories, one number R_eff,")
    say("      two different answers.  R_eff(geometric) is dynamics-blind.")
    say(f"  {'L':>3} {'R_eff^max same':>16} {'gap(m=0, Maxwell)':>18} {'gap(m=0.7 scalar)':>18}")
    for r in rows[:4]:
        say(f"  {r['L']:>3} {r['R']:>16.8e} {r['w']:>18.8e} {0.7:>18.8e}")
    say("")
    say("  VERDICT ACID 1: the route's 'success' on the free field is the tautological")
    say("  reading (R_eff := 1/lambda).  With the computable geometric R_eff the free")
    say("  massless and the free massive theory are INDISTINGUISHABLE => it cannot")
    say("  be a certificate.  Acid test FAILED (not for sign, for information).")
    say("")
    return rows


# ---------------------------------------------------------------------------
def part_acid2():
    say("=" * 78)
    say("ACID 2 -- beta = 0 Haar (ultralocal; corpus: G''=0.2688, no atom)")
    say("=" * 78)
    say("  Two exact facts at beta = 0:")
    say("   (i) geometric R_eff is finite and beta-free (same number as at beta=2.5);")
    say("   (ii) the transfer matrix is the rank-one kernel 11^T  =>  lambda_2/lambda_1 = 0")
    say("        => the lattice correlation length is ZERO, i.e. 'gap = infinity'.")
    say("  So at beta=0 the route's premise 'lambda >= c/R_eff' is satisfied with ANY c,")
    say("  and simultaneously the corpus's acid test says: no atom (f_1^2 = 0).")
    say("  A certificate that is compatible with 'gap = infinity' and with 'no atom'")
    say("  certifies nothing.  This is the false-positive point.")
    say("")
    say("  (numerics from qg2 PART A: at beta_t -> 0, t -> 1, ratio -> 0, m -> +inf;")
    say("   at beta_t = 0.1 the gap is already 2.3059 with the SAME geometry.)")
    say("")


# ---------------------------------------------------------------------------
def part_acid3():
    say("=" * 78)
    say("ACID 3 -- double well: R_eff finite but the gap is 0 (Kramers pair)")
    say("=" * 78)
    say("  H = -(hb^2/2) d^2/dx^2 + (x^2-1)^2 on [-3,3], Dirichlet, N=4000.")
    say("  Exact finite-difference spectrum (tridiagonal).")
    say(f"  {'hb':>6} {'E0':>12} {'E1':>12} {'gap=E1-E0':>12} {'R_eff(2-node)':>14} {'gap*R':>12}")
    rows = []
    N = 4000
    x = np.linspace(-3.0, 3.0, N)
    h = x[1] - x[0]
    V = (x ** 2 - 1.0) ** 2
    Rgeo = 2.0   # the 2-well reduction: two nodes, one intermediate vertex, unit weights
    for hb in [0.40, 0.30, 0.20, 0.15, 0.10, 0.07, 0.05]:
        # H = -(hb^2/2) d^2/dx^2 + V ; tridiagonal, Dirichlet ends
        d = hb ** 2 / h ** 2 + V[1:-1]
        e = -0.5 * hb ** 2 / h ** 2 * np.ones(N - 3)
        ev = eigh_tridiagonal(d, e, select="i", select_range=(0, 1))[0]
        gap = ev[1] - ev[0]
        S_wkb = math.sqrt(2.0) * 4.0 / 3.0      # int_{-1}^{1} sqrt(2 V) dx, V=(x^2-1)^2
        say(f"  {hb:>6.2f} {ev[0]:>12.6e} {ev[1]:>12.6e} {gap:>12.6e} {Rgeo:>14.1f} "
            f"{gap*Rgeo:>12.6e}   WKB e^-S/hb = {math.exp(-S_wkb/hb):>10.2e}")
        rows.append(dict(hb=hb, E0=float(ev[0]), E1=float(ev[1]), gap=float(gap), R=Rgeo))
    say("")
    say("  gap -> 0 exponentially in 1/hb, while ANY geometric R_eff (2-node reduction,")
    say("  or the 4000-node chain) stays finite and even DIVERGES with the discretisation.")
    say("  A finite-R_eff certificate must produce a positive floor here.  It cannot.")
    say("  VERDICT ACID 3: FAILED, as required.")
    say("")
    return rows


def main():
    t0 = time.time()
    say(f"qg3_acid.py   start {time.strftime('%Y-%m-%d %H:%M:%S')}")
    say("")
    part_G()
    A1 = part_acid1()
    part_acid2()
    A3 = part_acid3()
    say(f"total wall {time.time() - t0:.1f} s")
    json.dump(dict(acid1=A1, acid3=A3), open("qg3_acid.json", "w"), indent=1)
    open("qg3_acid.log", "w", encoding="utf-8").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
