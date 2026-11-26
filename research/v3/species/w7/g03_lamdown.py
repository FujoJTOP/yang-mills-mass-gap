"""G03 - Lambda_star's LOWER side: is "a lower bound for Lambda_star" the SAME proposition as I2?

    Lambda_star := inf supp( nu - f1^2 delta_{m1^2} ),   m1^2 = inf supp nu
    I2          : Lambda_star - m1^2 >= delta > 0        (uniform in (L,a,beta))

Part A  the only unconditional lower bound is the TRIVIAL one, Lambda_star >= m1^2 (vacuous)
Part B  thin gap (W4): identical (M0,M1,M2)=(2,11,101), isolation -> 0  => moments cannot lower-bound
Part C  Lemma C sharpness: Lambda_max = Lambda_star (2 atoms) / over-shoots (continuum); Lambda_max >= m1^2 always
Part D  the lower side IS attainable once the state count K is known (Prony exact) and in 2D YM (margin 5/4)
Part E  (Lambda_star >= m1^2+delta) <=> (no support in (m1^2, m1^2+delta))  -- the definitional identity = I2
Part F  a support point gives an UPPER bound: Lambda_star <= s_th (two-particle threshold)

    py -3.12 g03_lamdown.py > g03_lamdown.log
"""
import numpy as np

rng = np.random.default_rng(20261125)


def moments(xs, ws, k):
    xs, ws = np.asarray(xs, float), np.asarray(ws, float)
    return np.array([float(np.sum(ws * xs ** j)) for j in range(k)])


def lam_star(xs, ws, m2):
    """inf supp(nu - f1^2 delta_{m2}); +inf if the measure is a single atom at m2."""
    xs, ws = np.asarray(xs, float), np.asarray(ws, float)
    rest = xs[ws > 0]
    rest = rest[np.abs(rest - m2) > 1e-12]
    return float(rest.min()) if rest.size else np.inf


# ------------------------------------------------------------------ Part A
def part_A():
    print("== A: the sole unconditional lower bound is the trivial Lambda_star >= m1^2 ==")
    bad, gaps = 0, []
    for _ in range(20000):
        xs = np.sort(1.0 + rng.random(4) * 10)          # m1^2 = xs[0]
        ws = rng.random(4) + 1e-3
        gap = lam_star(xs, ws, xs[0]) - xs[0]
        bad += gap < -1e-12
        gaps.append(gap)
    gaps = np.array(gaps)
    print(f"   20000 random 4-atom measures: violations of Lambda_star >= m1^2 = {bad}")
    print(f"   gap: min {gaps.min():.3e}  median {np.median(gaps):.3f}  max {gaps.max():.3f}")
    print("   => [T] Lambda_star >= m1^2 always (0 violations); the sample min gap is >0 only because the")
    print("      support is distinct -- the exact infimum 0 is realised by the tower (Part B / Lemma B).")
    print("   => the trivial bound is all positivity alone gives; anything above it is an ABSENCE")
    print("      statement (Part E), not a bound of the same kind.\n")


# ------------------------------------------------------------------ Part B
def w4_row(eps, q_frac=0.4):
    """3 atoms {1, 1+eps, L} with q = q_frac*eps, reproducing (M0,M1,M2)=(2,11,101) in closed form."""
    q = q_frac * eps
    Lp1 = (99.0 - q * (2 * eps + eps * eps)) / (9.0 - q * eps)   # = L + 1
    L = Lp1 - 1.0
    r = (9.0 - q * eps) / (Lp1 - 2.0)                             # r(L-1) = 9 - q*eps
    p = 2.0 - q - r
    return np.array([1.0, 1.0 + eps, L]), np.array([p, q, r])


def part_B():
    print("== B: thin gap (W4 reproduced in closed form) -- fixed (M0,M1,M2)=(2,11,101) of nu0=d_1+d_10 ==")
    print("   reference nu0 = delta_1 + delta_10 : isolation = 9 exact, Lambda_max = 10 exact")
    print("       eps        L        p@1        q@(1+eps)   r@L       |M1-11|   |M2-101|   isolation")
    rows = []
    for eps in [0.5, 0.05, 0.01, 1e-3, 1e-6, 1e-9]:
        xs, ws = w4_row(eps)
        m = moments(xs, ws, 3)
        rows.append((eps, xs, ws, m))
        print(f"  {eps:8.1e} {xs[2]:9.5f} {ws[0]:10.6f} {ws[1]:11.3e} {ws[2]:10.6f}"
              f" {abs(m[1]-11):.1e} {abs(m[2]-101):.1e} {eps:10.1e}")
    print("   (eps=0.5 row reproduces R2-ISOLATION W4 verbatim: p=.82149 q=.2 r=.97851 L=10.0955)")
    print("   => [T]+[E] the fibre of the moment map at (2,11,101) reaches isolation eps for EVERY eps>0:")
    print("      inf over the fibre = 0. Exact realisation of 0: the tower/mu touching m1^2 (Lemma B).")
    print("      => no functional of finitely many moments lower-bounds Lambda_star - m1^2. [NG]\n")
    return rows


# ------------------------------------------------------------------ Part C
def part_C():
    print("== C: Lemma C (Lambda_star <= Lambda_max) -- sharp for 2 atoms, over-shoots with a continuum ==")
    m2, f2, S = 1.0, 0.6, 7.0
    M0 = 1.0
    M1 = f2 * m2 + (1 - f2) * S
    print(f"   2 atoms {{1,{S}}}, f1^2={f2}: Lambda_star={S:.4f}  Lambda_max={(M1-f2*m2)/(M0-f2):.4f}"
          f"  (equal => sharp)")
    s = m2 + np.linspace(0.0, 40.0, 400001)[1:]
    dens = np.exp(-(s - m2))
    Z = np.trapezoid(dens, s)
    tail_w = 1 - f2
    M1c = f2 * m2 + tail_w * np.trapezoid(s * dens, s) / Z
    M2c = f2 * m2 ** 2 + tail_w * np.trapezoid(s * s * dens, s) / Z
    lam_max = (M1c - f2 * m2) / (M0 - f2)
    print(f"   + continuum e^(-(s-1)) on (1,inf), weight {tail_w:.1f}: Lambda_star={m2:.4f} (touches m1^2),"
          f" Lambda_max={lam_max:.4f} > m1^2")
    print(f"   => Lemma C bound admits isolation 0 on the SAME (M0,M1)=(1,{M1c:.4f}) datum;")
    print(f"   monotone check: M1 - m1^2*M0 = {M1c - 1.0:+.4f} >= 0 => Lambda_max >= m1^2 always.")
    print(f"   (2-atom datum M1={M1:.3f}: Lambda_max={(M1-f2*m2)/(M0-f2):.4f} = Lambda_star, sharp)")
    print("   => [T] the bounds allow Lambda_star - m1^2 in [0,inf), infimum 0 => no finite margin.\n")


# ------------------------------------------------------------------ Part D
def prony(M, K):
    H = np.array([[M[n + j] for j in range(K)] for n in range(K)])
    rhs = -np.array([M[K + n] for n in range(K)])
    try:
        a = np.linalg.solve(H, rhs)
    except np.linalg.LinAlgError:
        return None
    roots = np.roots(np.concatenate([[1.0], a[::-1]]))   # prod(z-x_i) = z^K + a_{K-1}z^{K-1}+...+a_0
    if np.max(np.abs(roots.imag)) > 1e-6:
        return None
    xs = np.sort(roots.real)
    try:
        ws = np.linalg.solve(np.vander(xs, N=K, increasing=True), np.array(M[:K]))
    except np.linalg.LinAlgError:
        return None
    return xs, ws


def part_D():
    print("== D: the lower side IS attainable when the state count K is known (so NOT structurally dead) ==")
    xs0, ws0 = np.array([1.0, 2.5, 7.0]), np.array([0.5, 0.3, 0.2])
    out = prony(moments(xs0, ws0, 6), 3)
    print(f"   3-atom tower, K=3 given: recovered {np.round(out[0],10)}  isolation={out[0][1]-out[0][0]:.12f}"
          f"  err={np.max(np.abs(out[0]-xs0)):.2e}   [T] (Prony)")
    xs8 = 1.0 + np.arange(8) * 0.7
    ws8 = np.array([.30, .25, .18, .12, .07, .04, .025, .015])
    ws8 = ws8 / ws8.sum()
    p3 = prony(moments(xs8, ws8, 6), 3)
    if p3 is None:
        print(f"   8-atom continuum probed at the WRONG K=3: Hankel system singular / complex roots"
              f" => no K=3 measure matches the 6 moments  (true isolation {xs8[1]-xs8[0]:.4f})")
    else:
        print(f"   8-atom continuum probed at the WRONG K=3: 'isolation' {p3[0][1]-p3[0][0]:.4f}"
              f" vs true {xs8[1]-xs8[0]:.4f}  => meaningless")
    C2 = {1: 0.75, 2: 2.0}                       # SU(2) Casimir: fundamental 3/4, adjoint 2
    print(f"   2D YM SU(2) Migdal tower: m1^2=C2(f)={C2[1]}, Lambda_star=C2(adj)={C2[2]}"
          f" => margin {C2[2]-C2[1]} = 5/4 exact, rho^2={C2[2]/C2[1]:.6f} = 8/3   [M]/[E]")
    print("   => a non-trivial lower bound Lambda_star > m1^2 is a THEOREM in solvable models (exact margin).")
    print("      What is missing in 4D is the state count / uniformity, not the definition. No structural")
    print("      no-go on the definition; the dead route is moment/positivity (Part B), not the object.\n")


# ------------------------------------------------------------------ Part E
def part_E(rows):
    print("== E: (Lambda_star >= m1^2+delta) <=> (no support in (m1^2, m1^2+delta)) -- definitional ==")
    ok = tot = 0
    for eps, xs, ws, m in rows:
        for delta in [eps * 0.5, eps, eps * 2.0]:
            lb = (eps >= delta)
            empty = not np.any((xs > 1.0 + 1e-12) & (xs < 1.0 + delta - 1e-12))
            ok += (lb == empty)
            tot += 1
    print(f"   checked {tot} (measure, delta) pairs: both formulations agree {ok}/{tot}")
    print("   => [T] 'a lower bound strictly above m1^2' and 'no spectrum in (m1^2, Lambda_star)' are the")
    print("      SAME formula; with delta fixed uniformly in (L,a,beta) it is I2 verbatim.\n")


# ------------------------------------------------------------------ Part F
def part_F():
    print("== F: an exhibited support point is an UPPER bound (threshold carries no lower-bound content) ==")
    pts, ws = np.array([1.0, 4.0, 4.3]), np.array([0.5, 0.3, 0.2])
    print(f"   nu=.5 d_1+.3 d_4+.2 d_4.3: m1^2=1, continuum edge s_th=4, Lambda_star={lam_star(pts,ws,1.0):.1f}"
          f"  => Lambda_star <= s_th")
    print("   => [T] '2m1 in supp nu' / 's_th in supp nu' yields an UPPER bound on Lambda_star, never a lower one.\n")


if __name__ == "__main__":
    part_A()
    rows = part_B()
    part_C()
    part_D()
    part_E(rows)
    part_F()
    print("VERDICT  (a) fires: 'Lambda_star lower bound' == I2 verbatim => ONE missing piece, not two;")
    print("         (b) no separate piece: the only extra content is uniformity of delta = I2's own quantifier;")
    print("         (c) no structural no-go: non-trivial lower bounds exist as theorems (2D margin 5/4; Prony")
    print("             given K). What dies is the moment/positivity route [NG] - a method, not the object.")
    print("No QED. Clay 0%.")
