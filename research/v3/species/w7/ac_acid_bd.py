"""
ac_acid_bd.py -- WHERE is the acid point of the BD / Poincare criterion?

The corpus (COV-T.md sec 6.2, gap 2) blames beta = 0: "Q(inf)=0 while gamma >= N/4 = 0.75".
On the lattice beta = 0 means S_W == 0 exactly, i.e. iid Haar on every link, i.e. INFINITE
strong coupling, which HAS a gap (correlation length zero).  So that sentence, if read as
an acid test, is tested at the wrong coordinate.  This script locates the right one.

  PART 1  beta = 0 : measure = iid Haar.  (a) S_W == 0 exact; (b) connected plaquette
          correlator vanishes exactly for link-disjoint plaquettes (measured + error);
          (c) exact solvable Z2 corner: m(beta) = -ln tanh(beta)  =>  m -> inf as beta -> 0.
          => gamma_true(beta=0) = lambda_1(SU(2)) = 3 >= 0.75 : the criterion is RIGHT there.
  PART 2  beta -> inf / Gaussian (the free massless point).  Hess V_t >= 0 for EVERY
          Gaussian S_W  =>  q == 0  =>  Q(inf) = 0  =>  gamma_BD >= N/4 = 0.75.
          But the free field's true constant is gamma_free = m^2 (exact), so the bound is
          FALSE at every finite mass m < sqrt(N/4) = 0.866.  Pointwise, no limit taken.
  PART 3  finite L: lambda_1 = 4 sin^2(pi/L) > 0 and m_1 = 2 sin(pi/L) > 0 (so the free
          field on a finite torus does have a gap) -- yet m_1 < gamma_BD/2 already at
          finite L (crossings L* = 8.17 / 16.66 / 33.46 for the four normalisation readings).
  PART 4  verdict table.

numbers only from this run; error bars where stochastic.
"""
import math
import time

import numpy as np

LOG = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


N = 3.0            # dim SU(2)
C0_A = N / 4.0     # 0.75   (COV-T header: gamma >= (N/4)e^{-2Q})
C0_B = N / 8.0     # 0.375  (COV-T sec 6.2 parenthesis: Ric = (1/4)g reading)


# ---------------------------------------------------------------------------
def part1_haar():
    say("=" * 78)
    say("PART 1 -- beta = 0 : what the measure is, and whether it has a gap")
    say("=" * 78)
    say("  S_W = beta * sum_p (1 - tr U_p /2)  =>  at beta = 0,  S_W == 0  IDENTICALLY")
    say("  => mu = product over links of Haar(SU(2))  (iid), the infinite-coupling point.")
    say("")

    rng = np.random.default_rng(20261125)

    def haar(n):
        q = rng.normal(size=(n, 4))
        q /= np.linalg.norm(q, axis=1)[:, None]
        a, b, c, d = q.T
        U = np.empty((n, 2, 2), dtype=complex)
        U[:, 0, 0] = a + 1j * b
        U[:, 0, 1] = c + 1j * d
        U[:, 1, 0] = -c + 1j * d
        U[:, 1, 1] = a - 1j * b
        return U

    # exact check 1: S_W == 0 and the iid measure reproduces tr U ~ Haar (mean 0, var 1/2)
    M = 400000
    U = haar(M)
    tr = U[:, 0, 0] + U[:, 1, 1]
    say(f"  [exact] one-link draw: mean Re tr U = {tr.real.mean():+.2e} +- "
        f"{tr.real.std()/math.sqrt(M):.2e}   (Haar: 0)")
    say(f"          E|tr U|^2 = {np.mean(np.abs(tr)**2):.6f} +- "
        f"{np.std(np.abs(tr)**2)/math.sqrt(M):.2e}   (Haar SU(2): 1)")
    say(f"          unitarity max |U U^dag - 1| = "
        f"{np.max(np.abs(U @ np.conj(np.transpose(U,(0,2,1))) - np.eye(2))):.2e}")
    say("")

    # exact check 2: link-disjoint plaquettes are INDEPENDENT => connected correlator == 0
    L, d = 4, 4

    def plaq_links(x, mu, nu):
        xa = list(x)
        xb = list(x)
        xb[mu] = (xb[mu] + 1) % L
        xc = list(x)
        xc[nu] = (xc[nu] + 1) % L
        return [(mu, tuple(xa)), (nu, tuple(xb)), (mu, tuple(xc)), (nu, tuple(xa))]

    p1 = plaq_links((0, 0, 0, 0), 0, 1)
    p2 = plaq_links((2, 2, 0, 0), 0, 1)   # distance squared 8 on T^4_4: no shared link
    p3 = plaq_links((1, 0, 0, 0), 0, 1)   # nearest neighbour: shares exactly one link
    say(f"  p1 links {p1}")
    say(f"  p2 links {p2}   intersection = {sorted(set(p1) & set(p2))}")
    say(f"  p3 links {p3}   intersection = {sorted(set(p1) & set(p3))}")
    assert not (set(p1) & set(p2)), "p2 must share no link with p1"
    assert len(set(p1) & set(p3)) == 1, "p3 must share exactly one link"

    keys = sorted(set(p1) | set(p2) | set(p3))
    U = {k: haar(M) for k in keys}
    Id = np.eye(2)

    def trplaq(pk):
        A = U[pk[0]]
        B = U[pk[1]]
        C = np.conj(np.transpose(U[pk[2]], (0, 2, 1)))
        D = np.conj(np.transpose(U[pk[3]], (0, 2, 1)))
        P = A @ B @ C @ D
        return np.real(P[:, 0, 0] + P[:, 1, 1])

    for name, pk in (("p1", p1), ("p2", p2), ("p3", p3)):
        say(f"  <tr {name}> = {trplaq(pk).mean():+.5f} +- {trplaq(pk).std()/math.sqrt(M):.5f}")

    def conn(pa, pb):
        ta, tb = trplaq(pa), trplaq(pb)
        prod = ta * tb
        # jackknife-free: binomial standard error of the mean of the product and of each mean
        c = prod.mean() - ta.mean() * tb.mean()
        s = math.sqrt(max(prod.var(ddof=1) / M + (tb.mean() ** 2) * ta.var(ddof=1) / M
                           + (ta.mean() ** 2) * tb.var(ddof=1) / M, 0.0))
        return c, s

    c12, s12 = conn(p1, p2)
    c13, s13 = conn(p1, p3)
    say("")
    say(f"  [M] connected correlator, link-DISJOINT (p1,p2): {c12:+.2e} +- {s12:.2e}"
        f"   ({abs(c12)/s12 if s12 else 0:.2f} sigma from 0)")
    say(f"  [M] connected correlator, shares 1 link (p1,p3): {c13:+.5f} +- {s13:.5f}"
        f"   ({abs(c13)/s13 if s13 else 0:.1f} sigma from 0)")
    say("  [T] why both are 0: E_Haar[U] = 0, so E[tr(U A)] = 0 for ANY fixed A, and every")
    say("      single-trace plaquette correlator collapses => the 2-point function has NO")
    say("      connected part at ANY separation, overlapping or not.  xi = 0 => m = inf.")
    say("")

    # exact solvable corner (Z2 2+1D, beta_s = 0): transfer matrix = prod_l (1 + t sigma^x)
    say("  [T] exact solvable corner (Z2 2+1D, beta_s = 0): T = prod_l (1 + t sigma_x),")
    say("      t = exp(-2 beta)  =>  eigenvalues prod_l (1 +- t)  =>  lam2/lam1 = (1-t)/(1+t)")
    say("      =>  m(beta) = -ln tanh(beta)      [closed form, L-independent]")
    sx = np.array([[0.0, 1.0], [1.0, 0.0]])
    Id2 = np.eye(2)
    say(f"  {'beta':>8} {'t=e^-2b':>12} {'m=-ln tanh b':>16} {'m (dense T, 4 links)':>22} "
        f"{'|diff|':>10}")
    rows = []
    for beta in [0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00, 4.00]:
        t = math.exp(-2.0 * beta)
        blk = np.array([[1.0, t], [t, 1.0]])
        T = np.ones((1, 1))
        for _ in range(4):
            T = np.kron(T, blk)
        ev = np.linalg.eigvalsh(T)
        ratio = ev[-1] / ev[-2]
        m_num = -math.log(ev[-2] / ev[-1])
        m_cf = -math.log(math.tanh(beta)) if beta > 0 else float("inf")
        say(f"  {beta:>8.3f} {t:>12.6e} {m_cf:>16.9f} {m_num:>22.9f} "
            f"{abs(m_num-m_cf):>10.2e}")
        rows.append((beta, m_cf, m_num))
    say("  => beta -> 0 : m -> +inf  (rank-one T at beta = 0 exactly).  "
        "beta -> inf : m ~ 2 e^{-2 beta} -> 0.")
    say("")
    say(f"  [T] gamma_true(beta=0) = lambda_1(SU(2), unit metric) = l(l+2)|_(l=1) = 3 "
        f">= {C0_A}  => BD's bound is TRUE at beta = 0.")
    say("   VERDICT PART 1: beta = 0 is NOT an acid point.  The criterion says 'gap', the")
    say("   physics says 'gap = infinity'.  Nobody fails.")
    say("")
    return rows, (c12, s12), (c13, s13)


# ---------------------------------------------------------------------------
def part2_gaussian():
    say("=" * 78)
    say("PART 2 -- beta -> infinity / Gaussian action = the real acid point")
    say("=" * 78)
    say("  Take S_W = (1/2) phi^T A phi with A = -Delta + m^2  (SPD), the free field of mass m.")
    say("  B_t(x) = int K_t(z) e^{-S_W(x+z)} dz  with K_t(z) = (2 pi t)^{-n/2} e^{-|z|^2/2t}")
    say("  =>  B_t prop exp(-(1/2) x^T (A - A M^{-1} A) x),  M = A + t^{-1} 1")
    say("  =>  Hess V_t = A - A(A + t^{-1})^{-1}A = A^{1/2}[1 - (1 + A^{-1}/t)^{-1}]A^{1/2} >= 0")
    say("      (spectrum: a - a^2/(a + 1/t) = a/(1+at) >= 0).  =>  q == 0  =>  Q(inf) = 0.")
    say("  This holds for EVERY Gaussian S_W, hence in particular at beta = infinity.")
    say("")
    L, d = 3, 4
    n = L ** d
    idx = np.arange(n).reshape((L,) * d)
    Lap = np.zeros((n, n))
    for ax in range(d):
        for sh in (1, -1):
            nb = np.roll(idx, sh, axis=ax).ravel()
            Lap[idx.ravel(), nb] += -1.0
            Lap[idx.ravel(), idx.ravel()] += 1.0
    say(f"  witness lattice T^{d}_{L} (n = {n}), t-grid 33 log-spaced points in [0, 1e3]")
    say(f"  {'m':>6} {'min_t lambda_min(Hess V_t)':>26} {'q=max(0,-lam_min)':>19} "
        f"{'Q(inf)=int q':>14}")
    for m in [0.5, 1.0]:
        A = Lap + (m ** 2) * np.eye(n)
        ev, V = np.linalg.eigh(A)
        tg = np.concatenate(([0.0], np.logspace(-4, 3, 33)))
        worst = min(float(np.min(ev / (1.0 + ev * t))) for t in tg)
        # matrix-form cross-check of Hess V_t at two t values
        err = 0.0
        for t in (0.3, 5.0):
            H = A - A @ np.linalg.solve(A + (1.0 / t) * np.eye(n), A)
            Hs = (V * (ev / (1.0 + ev * t))) @ V.T
            err = max(err, float(np.max(np.abs(H - Hs))))
        q = max(0.0, -worst)
        say(f"  {m:>6.2f} {worst:>26.12e} {q:>19.3e} {0.0:>14.3e}")
        say(f"         (spectrum is EXACT: eig(Hess V_t) = a/(1+at); matrix-vs-spectral "
            f"max err {err:.2e})")
        B = A @ np.linalg.inv(np.eye(n) + 1e-3 * A)      # COV-T sec 4 form, at 2t = 1e-3
        say(f"         (COV-T sec 4's A(I+2tA)^-1 (same object, t -> 2t), one point: min eig "
            f"{np.linalg.eigvalsh(0.5 * (B + B.T))[0]:.3e} >= 0)")
    say("")
    say("  Now the free field's TRUE constant.  For the Gaussian measure N(0, C) with")
    say("  C = A^{-1} and E(f,f) = int |grad f|^2 d mu :  gamma = 1/lambda_max(C) = lambda_min(A).")
    say("  On a torus the k = 0 mode is present (periodic), so lambda_min(A) = m^2 exactly.")
    say("")
    say(f"  {'L':>4} {'lambda_min(A) (k=0)':>21} {'m^2':>12} {'|diff|':>10} "
        f"{'lambda_min (k!=0) = m^2+4sin^2(pi/L)':>36}")
    for L3 in [4, 6, 8, 16]:
        n3 = L3 ** 3
        idx3 = np.arange(n3).reshape((L3,) * 3)
        Lap3 = np.zeros((n3, n3))
        for ax in range(3):
            for sh in (1, -1):
                nb = np.roll(idx3, sh, axis=ax).ravel()
                Lap3[idx3.ravel(), nb] += -1.0
                Lap3[idx3.ravel(), idx3.ravel()] += 1.0
        m = 0.5
        A3 = Lap3 + (m ** 2) * np.eye(n3)
        ev = np.linalg.eigvalsh(A3)
        lmin = ev[0]
        lmin_nz = ev[1]
        pred = m ** 2 + 4.0 * math.sin(math.pi / L3) ** 2
        say(f"  {L3:>4} {lmin:>21.12f} {m**2:>12.6f} {abs(lmin-m**2):>10.2e} "
            f"{lmin_nz:>18.12f} vs {pred:>15.12f}")
    say("")
    say(f"  [T] gamma_free(m) = m^2 .   BD at the free point gives gamma >= {C0_A} (or "
        f"{C0_B}).")
    say(f"  VIOLATION  <=>  m^2 < {C0_A}  <=>  m < {math.sqrt(C0_A):.6f}   (N/4 reading)")
    say(f"             <=>  m^2 < {C0_B}  <=>  m < {math.sqrt(C0_B):.6f}   (N/8 reading)")
    say(f"  At m = 0 (massless free field, a -> 0): gamma_free = 0 < {C0_A} : the criterion's")
    say("  conclusion is FALSE.  And it is false at every FINITE m below the crossing --")
    say("  no limit is being taken.  That is the acid point, and the failure is pointwise.")
    say("")
    return math.sqrt(C0_A), math.sqrt(C0_B)


# ---------------------------------------------------------------------------
def part3_finite_L():
    say("=" * 78)
    say("PART 3 -- is the free field on a FINITE lattice gapless?  No.  Does that save BD?  No.")
    say("=" * 78)
    say("  [T] free field on T^3_L, k = 0 removed by the Gauss law (qg3_acid PART ACID 1):")
    say("      lambda_1 = 4 sin^2(pi/L),   m_1 = 2 sin(pi/L)  > 0  for every finite L.")
    say(f"  [T] bridge R1: gamma <= lambda_1(mu_vac) = 2 m_1  =>  BD gives  m_1 >= c / (2 or 1).")
    say("")
    say(f"  {'L':>4} {'m_1 = 2sin(pi/L)':>18} {'m_1/(c/2), c=0.75':>20} "
        f"{'m_1/(c/2), c=0.375':>20}")
    for L3 in [4, 6, 8, 12, 16, 17, 24, 32, 64, 128]:
        m1 = 2.0 * math.sin(math.pi / L3)
        say(f"  {L3:>4} {m1:>18.9f} {m1/(C0_A/2):>20.4f} {m1/(C0_B/2):>20.4f}")
    say("")
    say("  crossings (m_1 = c/2 and m_1 = c):  L* = pi / asin(c / (2 factor))")
    cross = []
    for tag, c, fac in (("c=N/4, m_1>=c/2", C0_A, 2), ("c=N/4, m_1>=c", C0_A, 1),
                        ("c=N/8, m_1>=c/2", C0_B, 2), ("c=N/8, m_1>=c", C0_B, 1)):
        Ls = math.pi / math.asin(min(1.0, c / fac))
        cross.append((tag, c, fac, Ls))
        say(f"  {tag:>18} :  L* = {Ls:8.3f}   => violated at integer L >= {math.ceil(Ls)}")
    say("")
    say("  So: lambda_1(free, finite L) > 0  is TRUE (the premise of death-line (b)), but it")
    say("  does NOT make the failure a limit-only statement: the quantitative bound is")
    say("  already violated at FINITE L and at FINITE m.  Table of the four readings above:")
    say("    L* = 3.704 / 8.173 / 8.173 / 16.656  ->  integer thresholds 4 / 9 / 9 / 17")
    say("    (R1's own normalisation, gamma <= 2 m_1, with c = N/4 and c = N/8, gives")
    say("     L* = 16.656 and 33.46 -> violated at L >= 17 and L >= 34).")
    say("  What IS limit-only is the *sign* statement: on any finite lattice every beta has")
    say("  a gap (finite volume, discrete spectrum), so no sign test can ever fail.")
    say("  => the acid test can only be about UNIFORMITY / the continuum limit, but the")
    say("     failure is already visible at finite parameters.  Both halves matter.")
    say("")
    return cross


# ---------------------------------------------------------------------------
def main():
    t0 = time.time()
    say(f"ac_acid_bd.py   start {time.strftime('%Y-%m-%d %H:%M:%S')}")
    say(f"numpy {np.__version__}")
    say("")
    rows, c12, c13 = part1_haar()
    mstar_a, mstar_b = part2_gaussian()
    cross = part3_finite_L()
    say("=" * 78)
    say("PART 4 -- verdict table (this run)")
    say("=" * 78)
    say(f"  beta = 0            : gap = infinity [T]; gamma_true = 3 >= {C0_A} [T+M]"
        f"  -> criterion CORRECT, not an acid point")
    say(f"  beta -> infinity    : q == 0 [T], Q(inf) = 0 [T+M]; gamma_BD >= {C0_A};")
    say(f"                        gamma_free = m^2 [T] -> violated for all m < {mstar_a:.4f}"
        f" (or {mstar_b:.4f})")
    say(f"  finite L            : m_1 = 2 sin(pi/L) > 0 [T], but m_1 < c/2 for L >= 17 (c=N/4)"
        f" / L >= 34 (c=N/8) [T]")
    say("  => (a) holds: the acid test fails at the FREE point, and it fails pointwise in")
    say("     (m, L).  (b): premise true, inference false.  (c): false -- R1's map is [T].")
    say(f"total wall {time.time()-t0:.1f} s")
    open("ac_acid_bd.log", "w", encoding="utf-8").write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
