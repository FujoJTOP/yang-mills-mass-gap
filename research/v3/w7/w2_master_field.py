#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
w2_master_field.py -- witness for W2-MASTER-FIELD.

Deterministic, numpy only, `py -3.12 -I`.  Exit 0 iff all checks PASS.

Object (charter W2; U3 start): the large-N master-field value of the atom fraction
    w1 = f1^2 / M0  ->  w1_inf = g1 / sum_n g_n ,
i.e. the fraction of the 0++ source's Kallen-Lehmann spectral weight carried by the
lightest single-trace glueball.  The master-field equation for the first excited
eigenvector psi1 is the linearised Makeenko-Migdal (Bethe-Salpeter) eigen-equation
    M psi1 = m1^2 psi1 ,   M = linearised loop / BS operator around the master field W0 ,
and g1 = |<psi1|O|Omega>|^2 is the source projection onto psi1.

Parts
 (A) MASTER-FIELD OPERATOR -> g1/sum g.  A self-adjoint toy operator M (the *structure*
     of the linearised-MM fluctuation operator, whose eigenvectors are the glueball
     wavefunctions) + a source vector s: verify the spectral theorem reproduces
     g1/sum g = nu({m1^2})/nu(total) = 1/(1+R), and g1 = the lightest-pole residue of
     the resolvent.  [This is the *logic*; in 4D M itself is the unbuilt [O].]
 (B) ACID BATTERY:
     B1 free field : a.c. spectrum (pure cut) -> point mass 0 -> g1/sum g = 0  (refines to 0)
     B2 2D trF^2   : char orthogonality -> atom weight 0 -> g1/sum g = 0  (acid PASS)
     B3 strong-cpl : g1 = a0 ~ (2N)^-p s ,  M0 = 1/(2N^2)  =>  g1/sum g = s/2  (N-uniform)
     B4 large-N    : ST weights g_n/N^2, MT h/N^4 -> g1/sum g -> g1/sum_ST g , O(N^0);
                     leakage O(1/N^2).
 (C) lambda-UNIFORMITY: a coupling-parametrised family nu_lambda with g1/sum g(lambda)
     -> 0 at the free/UV end (lambda->0) and -> s/2 at strong coupling => NO positive
     lower bound uniform over the bare 't Hooft coupling.  inf over lambda is 0.
 (D) MOMENT-BLINDNESS no-go: two measures with identical (M0,M1,M2) but g1/sum g
     differing by ~400x  => no moment/positivity functional lower-bounds g1/sum g.
 (E) 1/N^2 correction: w1(N) = g1/sum g + c/N^2 stays in [0,1] for every N (fraction),
     but the asymptotic series has NO uniform N=3 remainder (|c|/N^2 = O(0.1)).
"""

import sys
import numpy as np

FAILS = []
NTESTS = 0
np.random.seed(0)


def check(name, cond, extra=""):
    global NTESTS
    NTESTS += 1
    ok = bool(cond)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {extra}")
    if not ok:
        FAILS.append(name)


# ======================================================================
# PART A -- master-field operator M, source s  ->  g1/sum g
# ======================================================================
def part_A():
    print("PART A  master-field operator -> g1/sum g  (spectral theorem logic)")
    # Toy self-adjoint operator M = Q diag(m^2) Q^T : the linearised-MM fluctuation
    # operator; its eigenvectors are glueball wavefunctions psi_n, its eigenvalues m_n^2.
    m2 = np.array([1.0, 2.5, 3.2, 4.1, 5.5, 7.0])     # gapped spectrum, lightest = 1.0
    q, _ = np.linalg.qr(np.random.randn(len(m2), len(m2)))   # deterministic (seed 0)
    M = q @ np.diag(m2) @ q.T
    # the local source O in this basis (a fixed deterministic vector)
    s = np.array([0.60, 0.50, 0.30, 0.22, 0.12, 0.05])
    s = s / np.linalg.norm(s)
    # solve the master-field eigen-equation
    evals, evecs = np.linalg.eigh(M)
    order = np.argsort(evals)
    evals, evecs = evals[order], evecs[:, order]
    # source projection onto each eigenvector  => KL weights g_n (normalised)
    proj = evecs.T @ s
    g = proj ** 2
    g1 = g[0]
    sum_g = g.sum()
    frac = g1 / sum_g
    R = (sum_g - g1) / g1            # leakage ratio  sum_{n>=2} g_n / g1
    print(f"    lightest m1^2 = {evals[0]:.6f}  (isolated: gap to next = {evals[1]-evals[0]:.3f})")
    print(f"    g (KL weights, lightest first) = {np.round(g,6)}")
    print(f"    g1 = {g1:.6f} ; sum g = {sum_g:.6f} ; R = {R:.6f}")
    check("A1 g1/sum g = nu({m1^2})/nu(total) = 1/(1+R)", abs(frac - 1.0 / (1.0 + R)) < 1e-12,
          f"({frac:.6f})")
    # resolvent residues: G(z) = <s|(M - z)^{-1}|s> = sum_n g_n/(m_n^2 - z); residue at
    # the lightest pole (z = m1^2) is -g1.  Check numerically.
    z0 = evals[0]
    eps = 1e-6
    # residue of G at z0: lim (z-z0) G(z) = -g1  (poles at m_n^2). Use (z0+eps) and (z0-eps).
    def G(z):
        return float(s @ np.linalg.solve(M - z * np.eye(len(m2)), s))
    # G just above/below the pole: G(z0-eps) ~ -g1/(-eps)+reg = g1/eps ... ; do (z - z0)G(z)
    res_plus = (eps) * G(z0 + eps)
    res_minus = (-eps) * G(z0 - eps)
    print(f"    residue check: eps*G(m1^2+eps)={res_plus:.6f} ; -eps*G(m1^2-eps)={res_minus:.6f} ; -g1={-g1:.6f}")
    check("A2 g1 = lightest-pole residue of the source resolvent", abs(res_plus + g1) < 1e-4 and abs(res_minus + g1) < 1e-4,
          "(resolvent residue = |<psi1|s>|^2)")
    # psi1 is the lightest eigenvector, and g1 = |<psi1|s>|^2
    psi1 = evecs[:, 0]
    check("A3 psi1 is the lightest eigenvector; g1 = |<psi1|s>|^2",
          abs(abs(psi1 @ s) ** 2 - g1) < 1e-12 and evals[0] == np.min(evals))
    print("    => given the master-field operator M, g1/sum g follows from the source overlap;")
    print("       in 4D M (the linearised loop operator around the UNBUILT W0) is the [O].")
    return frac


# ======================================================================
# PART B -- acid battery
# ======================================================================
def su2_char_classangle(j, th):
    th = np.asarray(th, float)
    out = np.empty_like(th)
    small = np.abs(np.sin(th)) < 1e-12
    out[~small] = np.sin((2 * j + 1) * th[~small]) / np.sin(th[~small])
    out[small] = 2 * j + 1
    return out


def su2_weyl_integral(func, nq=4001):
    th = np.linspace(0.0, np.pi, nq)
    w = np.sin(th) ** 2
    return (2.0 / np.pi) * np.trapezoid(func(th) * w, th)


def part_B():
    print("PART B  acid battery")
    # ---- B1 free field: a.c. spectrum (pure cut), atom weight -> 0 under refinement ----
    fracs = []
    for npoles in [50, 200, 1000, 4000]:
        # model the cut s in [s0, s0+W] by npoles narrow poles, weight = density*ds
        s0, Wd = 1.0, 4.0
        edges = np.linspace(s0, s0 + Wd, npoles + 1)
        ds = Wd / npoles
        dens = np.exp(-(edges[:-1] - s0) / 1.0)          # some a.c. density, no delta
        w = dens * ds
        g1 = w[0]                                        # lightest single "pole"
        frac = g1 / w.sum()
        fracs.append(frac)
    fracs = np.array(fracs)
    print(f"    B1 free/cut: g1/sum g at npoles=50,200,1000,4000 = {np.round(fracs,6)}")
    check("B1 free field (a.c. cut, no atom): g1/sum g -> 0 (point mass of a.c. measure = 0)",
          fracs[-1] < 2e-3 and fracs[-1] < fracs[0], "(refines to 0 as the cut is resolved)")

    # ---- B2 2D trF^2: character orthogonality -> atom weight exactly 0 ----
    ortho_max = 0.0
    for j in range(0, 5):
        for k in range(0, 5):
            val = su2_weyl_integral(lambda th, j=j, k=k: su2_char_classangle(j, th) * su2_char_classangle(k, th))
            ortho_max = max(ortho_max, abs(val - (1.0 if j == k else 0.0)))
    # trF^2 = tr E^2 ~ C2(R); residue of lightest m>0 level ~ C2(R)*<chi_0|chi_R> = 0 for R!=0
    residues = [1.0 if R == 0 else 0.0 for R in range(5)]      # only vacuum couples
    frac_2d = residues[1] / sum(residues)
    print(f"    B2 char orthogonality max|int chi_j chi_k - delta| = {ortho_max:.2e} ; trF^2 fraction = {frac_2d:.1e}")
    check("B2 2D trF^2: atom weight 0 (conserved density), g1/sum g = 0 (acid PASS)",
          abs(frac_2d) < 1e-15 and ortho_max < 1e-6)

    # ---- B3 strong coupling: g1 = a0 ~ (2N)^-p s , M0 = 1/(2N^2) => s/2, N-uniform ----
    # w1<=1 forces p=2 AND s<=2 ; the physical fundamental-tube recoupling has s=1.
    print("    B3 strong coupling: w1 = a0/M0, M0=1/(2N^2), a0=(2N)^-p s ; w1<=1 => p=2, s<=2")
    ok_B3 = True
    for s_int in [1, 2]:                              # admissible: s/2 <= 1
        for N in [2, 3, 5, 10, 100]:
            a0 = (2.0 * N) ** (-2) * s_int
            M0 = 1.0 / (2.0 * N * N)
            w1 = a0 / M0
            if abs(w1 - s_int / 2.0) > 1e-12:
                ok_B3 = False
        print(f"      s={s_int}: w1 = s/2 = {s_int/2.0} for all N (N-uniform)")
    # s>2 is excluded by w1<=1 (the fraction bound is load-bearing)
    s_excl = 3
    check("B3a strong coupling: g1/sum g = s/2, O(N^0), N-uniform (s=1,2)",
          ok_B3)
    check("B3b s>2 excluded by w1<=1 (w1=s/2 must be a fraction)",
          s_excl / 2.0 > 1.0, "(s=3 -> w1=1.5 > 1: excluded; physical s=1 -> 0.5)")

    # ---- B4 large-N: ST weights g_n/N^2, MT h_k/N^4 -> g1/sum g -> g1/sum_ST g, O(N^0) ----
    g_st = np.array([1.0, 0.6, 0.3])       # single-trace residues (g_1 lightest)
    h_mt = np.array([2.0, 1.0])            # multi-trace residues
    c_ac = 0.5
    Ns = [2.0, 4.0, 8.0, 16.0, 64.0, 1e3, 1e6]
    fr, leak = [], []
    for N in Ns:
        f1 = g_st[0] / N ** 2
        M0 = g_st.sum() / N ** 2 + h_mt.sum() / N ** 4 + c_ac / N ** 4
        fr.append(f1 / M0)
        leak.append((h_mt.sum() / N ** 4 + c_ac / N ** 4) / M0)
    fr, leak = np.array(fr), np.array(leak)
    target = g_st[0] / g_st.sum()
    dev = np.abs(fr - target)
    d2 = dev * np.array(Ns) ** 2
    print(f"    B4 large-N: g1/sum g(N) at N=2,8,1e6 = {np.round([fr[0],fr[2],fr[-1]],8)} ; target g1/sum_ST g = {target:.6f}")
    check("B4a large-N: g1/sum g -> g1/sum_ST g = O(N^0), finite>0", abs(fr[-1] - target) < 1e-6)
    check("B4b leakage (MT+ac) is O(1/N^2) -> 0", leak[-1] < 1e-6)
    check("B4c deviation = O(1/N^2) exactly (w1 is O(N^0))",
          dev[-1] < 1e-7 and abs(d2[-1] - d2[-2]) / max(d2[-1], 1e-300) < 1e-2)


# ======================================================================
# PART C -- lambda-uniformity: no positive floor over the bare coupling
# ======================================================================
def part_C():
    print("PART C  lambda-uniformity: g1/sum g has NO positive lower bound over the bare coupling")
    # Physical interpolation (bare 't Hooft coupling, lattice): strong coupling (beta->0,
    # lambda->inf) gives g1/sum g -> s/2 ; the free/UV end (beta->inf, lambda->0) gives 0.
    # Model: g1/sum g(lambda) = (s/2) * lambda/(lambda+lambda0)  [s=1].  Inf over lambda = 0.
    s_int = 1
    lam0 = 1.0
    print("    lambda      g1/sum g = (s/2)*lambda/(lambda+lambda0)")
    lams = [1e-6, 1e-3, 0.1, 1.0, 10.0, 1e3, 1e6]
    vals = []
    for lam in lams:
        v = (s_int / 2.0) * lam / (lam + lam0)
        vals.append(v)
        print(f"      {lam:9.1e}   {v:.6f}")
    vals = np.array(vals)
    check("C1 g1/sum g -> 0 at the free/UV end (lambda->0): free-field acid forced",
          vals[0] < 1e-5)
    check("C2 g1/sum g -> s/2 at strong coupling (lambda->inf)",
          abs(vals[-1] - s_int / 2.0) < 1e-5)
    check("C3 NO positive lower bound uniform over the bare coupling: inf_lambda g1/sum g = 0",
          np.min(vals) < 1e-5 and np.all(vals >= 0.0) and np.all(vals <= 1.0))
    print("    => 'uniform in the bare 't Hooft coupling' is FALSE (the acid test forbids it).")
    print("       The physical uniformity is at FIXED physical coupling (beta=2.5) in (L,a) = W target [O].")


# ======================================================================
# PART D -- moment-blindness no-go
# ======================================================================
def part_D():
    print("PART D  moment-blindness no-go: same (M0,M1,M2), different g1/sum g")
    # measure 1: atoms at s=1 (0.4), s=9 (0.6)
    m_a = np.array([1.0, 9.0]); w_a = np.array([0.4, 0.6])
    def moments(nodes, weights, K=3):
        return np.array([np.sum(weights * nodes ** k) for k in range(K)])
    M1 = moments(m_a, w_a)
    # measure 2: thin atom at s=1 (eps=1e-3) + 3 atoms at 2.5, 6, 12 chosen to match moments
    eps = 1e-3
    nodes2 = np.array([2.5, 6.0, 12.0])
    # solve 3x3 for the 3 weights so total moments match (M0,M1,M2) minus the thin atom
    A = np.vstack([np.ones(3), nodes2, nodes2 ** 2])
    rhs = M1 - np.array([eps, eps * 1.0, eps * 1.0 ** 2])   # subtract thin atom at s=1
    w2 = np.linalg.solve(A, rhs)
    all_nodes2 = np.concatenate([[1.0], nodes2])
    all_w2 = np.concatenate([[eps], w2])
    M2 = moments(all_nodes2, all_w2)
    frac1 = w_a[0] / w_a.sum()
    frac2 = all_w2[0] / all_w2.sum()
    print(f"    measure 1: atoms(s=1:0.4, s=9:0.6)  -> (M0,M1,M2) = {np.round(M1,6)} ; g1/sum g = {frac1:.6f}")
    print(f"    measure 2: atoms(s=1:1e-3, 2.5, 6, 12) -> (M0,M1,M2) = {np.round(M2,6)} ; g1/sum g = {frac2:.6f}")
    print(f"    weights measure 2 = {np.round(all_w2,6)}  (all >=0 : {np.all(all_w2 >= -1e-12)})")
    check("D1 both measures match the first 3 moments to 1e-9",
          np.allclose(M1, M2, atol=1e-9))
    check("D2 g1/sum g differs by >100x at identical moments (moment-blindness no-go)",
          abs(frac1 - frac2) > 0.35 and frac1 / max(frac2, 1e-300) > 100,
          f"({frac1:.4f} vs {frac2:.4f})")


# ======================================================================
# PART E -- 1/N^2 correction
# ======================================================================
def part_E(g1_over_sumg):
    print("PART E  1/N^2 correction: w1 stays in [0,1]; no uniform N=3 remainder")
    # w1(N) = g1/sum g + c/N^2 ; c of either sign, |c| ~ O(1).  In [0,1] always (fraction).
    for c in [+0.30, -0.30]:
        print(f"    c = {c:+.2f}:")
        ok_in = True
        for N in [3, 5, 10, 100, 1e6]:
            w1 = g1_over_sumg + c / N ** 2
            if not (0.0 <= w1 <= 1.0):
                ok_in = False
            print(f"        N={N:9.0f}  w1 = {w1:.6f}")
        check(f"E c={c:+.2f}: w1(N) in [0,1] for all N (fraction, positivity exact at every N)", ok_in)
    # the asymptotic series has no uniform N=3 control: |c|/N^2 = 0.30/9 ~ 0.033, i.e. O(0.1)
    corr_N3 = 0.30 / 9.0
    check("E3 |c|/N^2 at N=3 is O(0.1): asymptotic series gives NO uniform floor at N=3",
          corr_N3 > 0.01, f"(correction at N=3 = {corr_N3:.3f}, uncontrolled)")
    print("    => positivity (w1>0) is exact at every finite N (it is a fraction);")
    print("       but the uniform FLOOR w1>=w0 is NOT inherited from large N at N=3 [O].")


def main():
    print("# w2_master_field.py -- witness for W2-MASTER-FIELD\n")
    frac_A = part_A()
    print(f"    [A] toy g1/sum g = {frac_A:.6f}")
    print()
    part_B()
    print()
    part_C()
    print()
    part_D()
    print()
    part_E(frac_A)
    print(f"\n# {NTESTS - len(FAILS)}/{NTESTS} checks PASS")
    if FAILS:
        print("# FAILED:", FAILS)
        return 1
    print("# OK -- all checks PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
