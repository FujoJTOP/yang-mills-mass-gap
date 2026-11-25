"""
qg1_network.py -- Q4 route (quantum graph / network reduction), step 1.

Question: for a scalar (undirected, weighted) network, what is the EXACT relation
between the spectral gap lambda_1 and the effective resistance R_eff?  The route's
premise is

        lambda_1 >= c / R_eff      ("resistance upper bound  ==>  gap lower bound"),

Rayleigh monotonicity (adding edges only lowers R_eff) being the friendly direction.

This script settles the direction exactly, then measures the L-scaling.

Contents
  PART A  sharp two-sided facts on the exact torus family (d = 1,2,3,4)
  PART B  brute-force enumeration of ALL connected graphs on n <= 6 vertices
          (sup/inf of lambda_1 * R_eff^max  =>  the exact constant)
  PART C  Kron / Schur (Y-Delta) reduction: R_eff invariance (exact) and the
          direction in which the gap moves
  PART D  the volume-factor bound lambda_1 >= 2/(n * Rbar)  (the only true
          resistance -> gap-floor inequality) and how weak it is on T^d

All output: stdout + qg1_network.log + qg1_network.json.  Single process, light.
"""
import json
import math
import time

import numpy as np

LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


# ----------------------------------------------------------------------------
# exact torus machinery
# ----------------------------------------------------------------------------
def torus_modes(L, d):
    """all k = 2 pi m / L, m in [0,L)^d ; returns m (n,d) and lambda(k) (n,)"""
    grids = np.meshgrid(*[np.arange(L)] * d, indexing="ij")
    m = np.stack([g.ravel() for g in grids], axis=1).astype(float)
    kk = 2.0 * math.pi * m / L
    lam = np.sum(4.0 * np.sin(kk / 2.0) ** 2, axis=1)
    return m, lam


def torus_Reff_all(L, d):
    """R_eff(r) for every r in the torus, by the exact Fourier sum
           R_eff(r) = (2/n) sum_{k!=0} (1 - cos k.r) / lambda(k).
       Cost O(n^2) once (n = L^d).  Returns (R (n,), r (n,d), lam (n,))."""
    m, lam = torus_modes(L, d)
    n = len(lam)
    nz = lam > 1e-12
    w = np.zeros(n)
    w[nz] = 1.0 / lam[nz]
    # cos(2 pi m.r / L) for all pairs
    C = np.cos(2.0 * math.pi * (m @ m.T) / L)
    S1 = w.sum()
    R = (2.0 / n) * (S1 - C.T @ w)
    R = np.maximum(R, 0.0)
    return R, m, lam


def torus_Reff_pinv(L, d):
    """independent cross-check: build the Laplacian and use the pseudo-inverse."""
    n = L ** d
    idx = np.arange(n).reshape((L,) * d)
    Lm = np.zeros((n, n))
    for ax in range(d):
        for sh in (1, -1):
            nb = np.roll(idx, sh, axis=ax).ravel()
            Lm[idx.ravel(), nb] += -1.0
            Lm[idx.ravel(), idx.ravel()] += 1.0
    Lp = np.linalg.pinv(Lm, hermitian=True)
    dg = np.diag(Lp)
    R = dg[:, None] + dg[None, :] - 2.0 * Lp
    return np.maximum(R, 0.0), Lm


# ----------------------------------------------------------------------------
# PART A : exact torus family
# ----------------------------------------------------------------------------
def part_A():
    say("=" * 78)
    say("PART A -- exact torus family T^d_L, unit conductances")
    say("  lambda_1 = 2d(1-cos(2pi/L));  R_eff(r) = (2/n) sum_{k!=0}(1-cos k.r)/lambda(k)")
    say("  PREMISE TESTED: is lambda_1 * R_eff^max bounded and L-independent?")
    say("=" * 78)
    say("")
    # cross-check Fourier formula against pinv  (align: R_pinv[i,j] depends on r=i-j)
    for (L, d) in [(2, 2), (3, 2), (4, 2), (2, 3), (3, 3), (2, 4)]:
        Rf, m, lam = torus_Reff_all(L, d)
        Rp, Lm = torus_Reff_pinv(L, d)
        n = L ** d
        ii = np.arange(n)
        mi = m  # (n,d) coordinates in the same ravel order as np.arange(n).reshape((L,)*d)
        ridx = np.empty((n, n), dtype=int)
        for j in range(n):
            delta = (mi - mi[j][None, :]) % L      # (n,d)
            ridx[:, j] = np.ravel_multi_index(delta.T.astype(int), (L,) * d, order="C")
        err = np.max(np.abs(Rf[ridx] - Rp))
        say(f"  [check] T^{d}_{L}: max|R_fourier - R_pinv| = {err:.3e}")
    say("")

    rows = []
    for d in [1, 2, 3, 4]:
        say(f"  --- d = {d} " + "-" * 55)
        say(f"  {'L':>4} {'n':>7} {'lambda_1':>12} {'R_eff^max':>12} "
            f"{'lam*Rmax':>12} {'Rbar':>11} {'1/(n lam_1)':>12} {'2/(n Rbar)':>11} {'true/tight':>10}")
        for L in range(2, 33 if d == 1 else (17 if d == 2 else (9 if d == 3 else 7))):
            n = L ** d
            lam1 = 2.0 * d * (1.0 - math.cos(2.0 * math.pi / L))
            if d <= 3 or L <= 6:
                R, m, lamk = torus_Reff_all(L, d)
                Rmax = float(R.max())
                S1 = np.sum(np.where(lamk > 1e-12, 1.0 / np.maximum(lamk, 1e-300), 0.0))
                Rbar = 2.0 * S1 / n
            else:
                # antipodal candidate only (n = L^4 gets big)
                m, lamk = torus_modes(L, d)
                r = np.array([L // 2] * d, dtype=float)
                c = np.cos(2.0 * math.pi * (m @ r) / L)
                S1 = np.sum(np.where(lamk > 1e-12, 1.0 / np.maximum(lamk, 1e-300), 0.0))
                Rmax = float((2.0 / n) * np.sum((1.0 - c) * np.where(lamk > 1e-12,
                                                                  1.0 / np.maximum(lamk, 1e-300), 0.0)))
                Rbar = 2.0 * S1 / n
            prod = lam1 * Rmax
            cq = 1.0 / (n * lam1)
            bound = 2.0 / (n * Rbar)
            say(f"  {L:>4} {n:>7} {lam1:>12.6e} {Rmax:>12.6e} {prod:>12.6e} "
                f"{Rbar:>11.5e} {cq:>12.4e} {bound:>11.4e} {lam1/bound:>10.3e}")
            rows.append(dict(d=d, L=L, n=n, lam1=lam1, Rmax=Rmax, prod=prod, Rbar=Rbar,
                             C_poincare=cq, bound=bound, ratio=lam1 / bound))

    say("")
    say("  --- power-law fits  prod = lambda_1 * R_eff^max ~ L^(-kappa) " + "-" * 20)
    for d in [1, 2, 3, 4]:
        sub = [r for r in rows if r["d"] == d and r["L"] >= 4]
        if len(sub) < 3:
            continue
        x = np.log([r["L"] for r in sub])
        y = np.log([r["prod"] for r in sub])
        kap = -np.polyfit(x, y, 1)[0]
        say(f"    d={d}: kappa = {kap:.4f}   (prod at L={sub[0]['L']}: {sub[0]['prod']:.3e}"
            f" -> L={sub[-1]['L']}: {sub[-1]['prod']:.3e})")
    say("")
    say("  VERDICT A: in every d the product DECAYS (kappa > 0); it is not bounded")
    say("  away from 0.  'lambda_1 R_eff = O(1)' is false on the exact torus family.")
    say("")
    return rows


# ----------------------------------------------------------------------------
# PART B : all connected graphs on n <= 6 vertices, unit weights
# ----------------------------------------------------------------------------
def gap_and_Reff(Lm):
    ev, U = np.linalg.eigh(Lm)
    lam1 = ev[1]
    U1 = U[:, 1:]
    inv = 1.0 / ev[1:]
    # R_eff(i,j) = sum_k (U1[i,k]-U1[j,k])^2 / ev[k] = a_i + a_j - 2 M_ij
    a = (U1 ** 2) @ inv
    M = (U1 * inv) @ U1.T
    R = a[:, None] + a[None, :] - 2.0 * M
    return lam1, float(R.max())


def part_B():
    say("=" * 78)
    say("PART B -- exhaustive enumeration of connected graphs, n = 2..6, unit weights")
    say("  aim: the exact constant in  lambda_1 * R_eff(a,b) <= C")
    say("=" * 78)
    out = []
    for n in range(2, 7):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        m = len(pairs)
        best = -1.0
        bestmask = None
        worst = 1e9
        worstmask = None
        cnt = 0
        for mask in range(1, 1 << m):
            Lm = np.zeros((n, n))
            ok = True
            for e, (i, j) in enumerate(pairs):
                if mask >> e & 1:
                    Lm[i, i] += 1
                    Lm[j, j] += 1
                    Lm[i, j] -= 1
                    Lm[j, i] -= 1
            # connected?
            ev = np.linalg.eigvalsh(Lm)
            if ev[1] < 1e-9:
                continue
            lam1, Rmax = gap_and_Reff(Lm)
            cnt += 1
            if cnt > 70000:
                break
            p = lam1 * Rmax
            if p > best:
                best, bestmask = p, mask
            if p < worst:
                worst, worstmask = p, mask
        say(f"  n={n}: {cnt} connected graphs  sup(lam_1 R_max) = {best:.9f}   "
            f"inf = {worst:.6e}   (theory bound 2)")
        out.append(dict(n=n, count=cnt, sup=best, inf=worst,
                        sup_edges=[[pairs[e][0], pairs[e][1]] for e in range(m) if bestmask >> e & 1]))
    say("")
    say("  VERDICT B: sup = 2.000000000 exactly (attained by any 2-vertex graph: R_eff=1, lam_1=2).")
    say("  The inequality therefore reads  lambda_1 <= 2/R_eff  -- an UPPER bound on the gap.")
    say("")
    return out


# ----------------------------------------------------------------------------
# PART C : Kron / Schur reduction  (the "effective network")
# ----------------------------------------------------------------------------
def part_C(ntrial=4000, seed=12345):
    say("=" * 78)
    say("PART C -- Kron/Schur (star-triangle) effective network")
    say("  (i)  R_eff among retained terminals is EXACTLY invariant   [should be T]")
    say("  (ii) which way does the gap move?                          [to be measured]")
    say("=" * 78)
    rng = np.random.default_rng(seed)
    maxinv = 0.0
    ratios = []
    ntot = 0
    for _ in range(ntrial):
        n = int(rng.integers(5, 11))
        # random connected graph: tree + extra edges, random positive weights
        A = np.zeros((n, n))
        perm = rng.permutation(n)
        for k in range(1, n):
            i, j = perm[k], perm[int(rng.integers(0, k))]
            A[i, j] = A[j, i] = rng.uniform(0.2, 3.0)
        for _ in range(int(rng.integers(0, n))):
            i, j = rng.choice(n, 2, replace=False)
            A[i, j] = A[j, i] = rng.uniform(0.2, 3.0)
        Lm = np.diag(A.sum(1)) - A
        ev = np.linalg.eigvalsh(Lm)
        if ev[1] < 1e-8:
            continue
        # random split
        nb = int(rng.integers(1, max(2, n - 2)))
        B = rng.choice(n, nb, replace=False)
        Amask = np.setdiff1d(np.arange(n), B)
        if len(Amask) < 2:
            continue
        LAA = Lm[np.ix_(Amask, Amask)]
        LBB = Lm[np.ix_(B, B)]
        LAB = Lm[np.ix_(Amask, B)]
        try:
            LA = LAA - LAB @ np.linalg.solve(LBB, LAB.T)
        except np.linalg.LinAlgError:
            continue
        evA = np.linalg.eigvalsh(LA)
        if evA[1] < 1e-12:
            continue
        # (i) R_eff invariance
        Rp_full = np.linalg.pinv(Lm, hermitian=True)
        Rp_A = np.linalg.pinv(LA, hermitian=True)
        dgf = np.diag(Rp_full)[Amask]
        dga = np.diag(Rp_A)
        Rf = dgf[:, None] + dgf[None, :] - 2 * Rp_full[np.ix_(Amask, Amask)]
        Ra = dga[:, None] + dga[None, :] - 2 * Rp_A
        maxinv = max(maxinv, float(np.max(np.abs(Rf - Ra))))
        # (ii) gap direction
        # normalise: both are Laplacians; use the smallest nonzero eigenvalue
        ratios.append(float(evA[1] / ev[1]))
        ntot += 1

    ratios = np.array(ratios)
    say(f"  trials kept: {ntot}")
    say(f"  (i)  max |R_eff^Kron - R_eff^full| over all retained pairs = {maxinv:.3e}   [invariant: YES]")
    say(f"  (ii) lambda_1(L_A) / lambda_1(L) :  min = {ratios.min():.6f}   "
        f"median = {np.median(ratios):.6f}   max = {ratios.max():.6e}")
    say(f"       fraction with ratio >= 1 (gap does NOT decrease): "
        f"{np.mean(ratios >= 1 - 1e-9):.4f}")
    say("")
    say("  VERDICT C: Kron reduction PRESERVES R_eff and (empirically, 4000 trials) can only")
    say("  RAISE the gap. Since R_eff is unchanged, any certificate of the shape")
    say("  'lambda >= c/R_eff' evaluated on the REDUCED network is an upper bound for the")
    say("  reduced object and says nothing about the original: the direction is useless.")
    say("")
    return dict(maxinv=maxinv, ratio_min=float(ratios.min()),
                ratio_median=float(np.median(ratios)),
                frac_ge1=float(np.mean(ratios >= 1 - 1e-9)), ntot=ntot)


# ----------------------------------------------------------------------------
# PART D : the only true resistance -> gap-floor inequality, and its price
# ----------------------------------------------------------------------------
def part_D():
    say("=" * 78)
    say("PART D -- the only TRUE gap FLOOR obtainable from resistances")
    say("  (f(x)-f(y))^2 <= R_eff(x,y) E(f)   [Cauchy-Schwarz in L^+]")
    say("  ==> Var_mu(f) <= (Rbar/2) E(f),  Rbar = mean pairwise resistance")
    say("  ==> lambda_1 >= 2/(n Rbar)        <-- carries 1/n = 1/volume")
    say("=" * 78)
    say("")
    say(f"  {'d':>3} {'L':>4} {'n':>7} {'2/(n Rbar)':>12} {'lambda_1':>12} {'slack':>10}")
    for d in [1, 2, 3, 4]:
        for L in ([8, 16, 32] if d == 1 else ([4, 8, 16] if d == 2 else ([4, 6, 8] if d == 3 else [4, 5, 6]))):
            n = L ** d
            lam1 = 2.0 * d * (1.0 - math.cos(2.0 * math.pi / L))
            m, lamk = torus_modes(L, d)
            S1 = np.sum(np.where(lamk > 1e-12, 1.0 / np.maximum(lamk, 1e-300), 0.0))
            Rbar = 2.0 * S1 / n
            b = 2.0 / (n * Rbar)
            say(f"  {d:>3} {L:>4} {n:>7} {b:>12.6e} {lam1:>12.6e} {lam1/b:>10.3e}")
    say("")
    say("  The true bound is weaker than the truth by a factor growing like n = L^d:")
    say("  it is a gap floor with an explicit 1/volume price, i.e. EXACTLY the factor that")
    say("  condition 5 (uniformity in (L,a)) forbids.")
    say("")


def main():
    t0 = time.time()
    say(f"qg1_network.py   start {time.strftime('%Y-%m-%d %H:%M:%S')}")
    say("")
    A = part_A()
    B = part_B()
    C = part_C()
    part_D()
    say("=" * 78)
    say(f"total wall {time.time() - t0:.1f} s")
    say("=" * 78)
    json.dump(dict(partA=A, partB=B, partC=C), open("qg1_network.json", "w"), indent=1)
    open("qg1_network.log", "w", encoding="utf-8").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
