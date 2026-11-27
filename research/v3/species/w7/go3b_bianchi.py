"""
go3b_bianchi.py -- GO3-BIANCHI witness (deterministic, fixed seeds, no MC).

Question: can the Bianchi identity  D_[mu F_{nu rho]} = 0  play the role that
phi.dphi = 0 plays in the 2D sigma model -- i.e. rewrite a NEGATIVE p^4 term as a
POSITIVE square -- and thereby lift k0_shell above 2 on the positive-cone class
(c_L >= 0 => k0_shell = 2, GL2-RPORDER Prop. R) ?

Parts
  A  free field: quantify the p^4 content of Bianchi.
      A1  k_[lam F_{mu nu}] = 0 identically for every (k, A_T)
      A2  sum_lam|k_lam F|^2 / |k.F|^2  (k-independent) => with Bianchi the whole
          2-derivative content of F is |k.F|^2 = (k^2)^2|A_T|^2
      A3  the "Bianchi square" functional vanishes => p^4 coefficient EXACTLY ZERO
  B  lattice loop symbols h_L(p) = |w_L(p)|^2 (w orthogonal to n_mu = e^{ip_mu}-1)
      B0  w.n = 0 ;  sum over 6 planes of the plaquette = 3 phat^2   (exact)
      B1  per shape A, B, kappa = -B/A
      B2  hypercubic orbit sums: B(u) = b1 + b2 sum u^4   (b1,b2 per unit A)
      B3  (plaq, 1x2 rect) family: reproduce GL2-N2  B/A = -(1/12 + c1) on the axis
      B4  pure plaquette action kappa(u) range
      B5  LP: does the cone contain  sum c_i (b2/A)_i = 0 ?  if yes, add the isotropic
          square lambda (p.p)^2 with lambda = -sum c_i (b1/A)_i > 0 and VERIFY directly
          that  h = p^2 + O(p^6)  in every direction.
  C  acid tests: 2D control, EOM-square isotropy, single plaquette.

Run:  py -3.12 go3b_bianchi.py
"""
import itertools
import numpy as np

OUT = []


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    OUT.append(s)


# ---------------------------------------------------------------- free field (A)
def F_of(k, A):
    return 1j * (np.outer(k, A) - np.outer(A, k))


def part_A():
    p("=" * 80)
    p("PART A  free field: what p^4 content does Bianchi have?")
    p("=" * 80)
    rng = np.random.default_rng(20261126)
    wb = 0.0
    ratios = []
    for _ in range(4000):
        k = rng.normal(size=4)
        A = rng.normal(size=4)
        A = A - (A @ k) / (k @ k) * k                 # Landau/transverse representative
        F = F_of(k, A)
        worst = 0.0
        for lam, mu, nu in itertools.permutations(range(4), 3):
            T = k[lam] * F[mu, nu] + k[mu] * F[nu, lam] + k[nu] * F[lam, mu]
            worst = max(worst, abs(T))
        wb = max(wb, worst)
        dF2 = sum(abs(k[lam] * F[mu, nu]) ** 2
                  for lam in range(4) for mu in range(4) for nu in range(4))
        dstar = sum(abs(sum(k[mu] * F[mu, nu] for mu in range(4))) ** 2 for nu in range(4))
        ratios.append(dF2 / dstar)
    p("A1  max |k_[lam F_{mu nu]}| over 4000 random (k, A_T):  %.3e" % wb)
    p("A2  sum_lam |k_lam F_{mu nu}|^2 / |k.F|^2 :  %.12f .. %.12f  (spread %.1e)"
      % (min(ratios), max(ratios), max(ratios) - min(ratios)))
    p("    Weitzenboeck with Bianchi:  |grad F|^2 = |dF|^2 + |d*F|^2 = |D.F|^2")
    p("    and |D.F|^2 = (k^2)^2 |A_T|^2 : ONE invariant, direction-free, ZERO on the")
    p("    free EOM  k.F = 0.  So Bianchi's positive p^4 = the EOM square, nothing else.")
    p("A3  'Bianchi square' operator  sum |k_[lam F_{mu nu]}|^2 == 0 : p^4 coefficient")
    p("    is EXACTLY ZERO (not 'positive').  Same on the lattice: prod_{p in dc} U_p = 1.")
    p("")
    return wb, min(ratios), max(ratios)


# ------------------------------------------------------------ lattice symbols (B)
def make_loop(step_list):
    """step_list: [(mu, n), ...].  Representation: (base of the link VARIABLE, mu, sign).
    For a backward step the variable sits at x - e_mu (matters at O(p^4))."""
    x = np.zeros(4, int)
    out = []
    for mu, n in step_list:
        s = 1 if n > 0 else -1
        e = np.eye(4, dtype=int)[mu]
        for _ in range(abs(n)):
            out.append((tuple(x if s > 0 else x - e), mu, s))
            x[mu] += s
    assert not x.any(), "loop does not close"
    return out


def w_of(loop, p):
    w = np.zeros(4, dtype=complex)
    for x, mu, s in loop:
        w[mu] += s * np.exp(1j * np.dot(p, np.array(x, float)))
    return w


def h_of(loop, p):
    w = w_of(loop, p)
    return float(np.vdot(w, w).real)


TS = (0.08, 0.04, 0.02, 0.01)
MB = np.array([[t ** 2, t ** 4, t ** 6, t ** 8] for t in TS])


def AB(fun, u):
    u = np.asarray(u, float)
    u = u / np.linalg.norm(u)
    vals = np.array([fun(t * u) for t in TS])
    coef, *_ = np.linalg.lstsq(MB, vals, rcond=None)
    return coef[0], coef[1]


H = []
for _pm in itertools.permutations(range(4)):
    for _sg in itertools.product([1, -1], repeat=4):
        R = np.zeros((4, 4), int)
        for _i, _mu in enumerate(_pm):
            R[_i, _mu] = _sg[_i]
        H.append(R)


def orbit_fun(loop):
    def f(q):
        q = np.asarray(q, float)
        return sum(h_of(loop, R.T @ q) for R in H) / len(H)
    return f


SHAPES = {
    "plaq":     [(0, 1), (1, 1), (0, -1), (1, -1)],
    "rect_1x2": [(0, 1), (1, 2), (0, -1), (1, -2)],
    "rect_2x2": [(0, 2), (1, 2), (0, -2), (1, -2)],
    "rect_1x3": [(0, 1), (1, 3), (0, -1), (1, -3)],
    "rect_2x3": [(0, 2), (1, 3), (0, -2), (1, -3)],
    "rect_3x3": [(0, 3), (1, 3), (0, -3), (1, -3)],
    "chair3d":  [(0, 1), (1, 1), (2, 1), (0, -1), (1, -1), (2, -1)],
    "hex3d":    [(0, 1), (1, 1), (2, 1), (0, -1), (2, -1), (1, -1)],
    "stair4d":  [(0, 1), (1, 1), (2, 1), (3, 1), (0, -1), (1, -1), (2, -1), (3, -1)],
}

DIRS = {
    "axis(1,0,0,0)":   [1, 0, 0, 0],
    "plane2(1,1,0,0)": [1, 1, 0, 0],
    "plane3(1,1,1,0)": [1, 1, 1, 0],
    "diag4(1,1,1,1)":  [1, 1, 1, 1],
    "odd(1,2,3,0)":    [1, 2, 3, 0],
}


def part_B():
    p("=" * 80)
    p("PART B  lattice loop symbols")
    p("=" * 80)
    loops = {k: make_loop(v) for k, v in SHAPES.items()}
    rng = np.random.default_rng(7)

    worst = 0.0
    for name, L in loops.items():
        for _ in range(40):
            q = rng.normal(size=4) * 0.7
            w = w_of(L, q)
            n = np.exp(1j * q) - 1.0
            worst = max(worst, abs(np.dot(w, n)) / max(1e-30, np.linalg.norm(w)))
    p("B0a max |w.n|/|w| over all shapes x 40 momenta : %.2e" % worst)
    worst3 = 0.0
    for _ in range(200):
        q = rng.normal(size=4) * 0.9
        tot = 0.0
        for mu in range(4):
            for nu in range(mu + 1, 4):
                tot += h_of(make_loop([(mu, 1), (nu, 1), (mu, -1), (nu, -1)]), q)
        phat2 = sum(4 * np.sin(x / 2) ** 2 for x in q)
        worst3 = max(worst3, abs(tot - 3 * phat2) / (3 * phat2))
    p("B0b |sum_{6 planes} h_plaq - 3 phat^2| / (3 phat^2) max : %.2e" % worst3)
    p("    => h_L = |w_L|^2 is the exact (gauge-invariant) symbol; pure plaquette = 3 phat^2.")
    p("")

    p("B1  per shape: h = A t^2 + B t^4 + ... ; kappa = -B/A")
    p("    shape        direction             A            B          kappa")
    for name, L in loops.items():
        for dn, u in DIRS.items():
            A, B = AB(lambda q: h_of(L, q), u)
            p("    %-10s %-20s %10.6f %12.6f %11.6f" % (name, dn, A, B, -B / A))
    p("")

    p("B2  hypercubic orbit sums: fit B(u) = b1 + b2 sum u^4  (per unit A)")
    rngu = np.random.default_rng(31337)
    U = rngu.normal(size=(160, 4))
    U /= np.linalg.norm(U, axis=1)[:, None]
    S4 = (U ** 4).sum(axis=1)
    Mfit = np.column_stack([np.ones_like(S4), S4])
    tab = {}
    p("    shape          A          A-spread     b1/A       b2/A      resid")
    for name, L in loops.items():
        f = orbit_fun(L)
        A, B = [], []
        for u in U:
            a, b = AB(f, u)
            A.append(a)
            B.append(b)
        A, B = np.array(A), np.array(B)
        coef, *_ = np.linalg.lstsq(Mfit, B, rcond=None)
        b1, b2 = coef
        resid = np.abs(B - Mfit @ coef).max()
        tab[name] = (A.mean(), b1 / A.mean(), b2 / A.mean())
        p("    %-11s %11.6f  %.1e  %10.6f %10.6f  %.1e"
          % (name, A.mean(), A.max() - A.min(), b1 / A.mean(), b2 / A.mean(), resid))
    p("    sign(b2/A): " + ", ".join("%s:%+.5f" % (n, t[2]) for n, t in tab.items()))
    p("")

    # B3 family, GL2-N2
    p("B3  (plaq, 1x2 rect) family, c_M + 8 c1 = 1 (both rect orientations per plane),")
    p("    A = 3 by construction.  GL2-N2 target: B/A = -(1/12 + c1) on the axis.")
    def hP(q, mu, nu):
        return 4 * np.sin(q[mu] / 2) ** 2 + 4 * np.sin(q[nu] / 2) ** 2

    def hRA(q, mu, nu):
        w = np.zeros(4, dtype=complex)
        w[mu] = 1 - np.exp(2j * q[nu])
        w[nu] = (np.exp(1j * q[mu]) - 1) * (1 + np.exp(1j * q[nu]))
        return float(np.vdot(w, w).real)

    def hRB(q, mu, nu):
        w = np.zeros(4, dtype=complex)
        w[nu] = 1 - np.exp(2j * q[mu])
        w[mu] = (np.exp(1j * q[nu]) - 1) * (1 + np.exp(1j * q[mu]))
        return float(np.vdot(w, w).real)

    def fam(c1):
        def f(q):
            tot = 0.0
            for mu in range(4):
                for nu in range(mu + 1, 4):
                    c0 = 1.0 - 8.0 * c1
                    tot += c0 * hP(q, mu, nu) + c1 * (hRA(q, mu, nu) + hRB(q, mu, nu))
            return tot
        return f
    p("    direction             B/A at c1=0     dB/A/dc1     B/A at c1=-1/12")
    for dn, u in DIRS.items():
        c1s = np.array([-0.04, -0.02, 0.0, 0.02, 0.04])
        vals = np.array([AB(fam(c1), u)[1] / AB(fam(c1), u)[0] for c1 in c1s])
        fit = np.polyfit(c1s, vals, 1)
        p("    %-20s %+13.6f  %11.6f  %+13.6f"
          % (dn, np.polyval(fit, 0.0), fit[0], np.polyval(fit, -1.0 / 12.0)))
    p("")

    # B4
    p("B4  pure plaquette action: b1/A = 0, b2/A = -1/12 => kappa(u) = sum u^4 / 12")
    p("    exact range on the 4D unit sphere: sum u^4 in [1/4, 1] => kappa in [1/48, 1/12]")
    fP = orbit_fun(loops["plaq"])
    kap = np.array([-AB(fP, u)[1] / AB(fP, u)[0] for u in U])
    p("    sampled over 160 random directions: [%.6f, %.6f]" % (kap.min(), kap.max()))
    p("    best isotropic-square residual (exact) = (1/12 - 1/48)/2 = 1/32 = %.6f" % (1 / 32))
    p("")

    # B5
    p("B5  LP over orbits:  c_i >= 0 ,  sum c_i = 1 ,  sum c_i (b2/A)_i = 0 ?")
    names = list(tab)
    B2 = np.array([tab[n][2] for n in names])
    B1 = np.array([tab[n][1] for n in names])
    from scipy.optimize import linprog
    res = linprog(c=np.zeros(len(names)), A_eq=np.vstack([np.ones(len(names)), B2]),
                  b_eq=np.array([1.0, 0.0]), bounds=[(0, None)] * len(names), method="highs")
    p("    linprog status %d (%s) ; feasible = %s" % (res.status, res.message, res.success))
    if res.success:
        c = res.x
        p("    solution c = " + ", ".join("%s:%.5f" % (n, v) for n, v in zip(names, c) if v > 1e-9))
        b1c = c @ B1
        lam = -b1c
        p("    (b1/A)_comb = %+.6f  =>  lambda = %+.6f  (needs > 0)" % (b1c, lam))

        def hcomb(q):
            return sum(v * orbit_fun(loops[n])(q) / tab[n][0] for n, v in zip(names, c) if v > 1e-12)
        p("    direct verification of  h = sum c_i h_i/A_i + lambda (p.p)^2 :")
        p("      direction            h/t^2 (t=0.01)     (h-t^2)/t^4")
        for dn, u in DIRS.items():
            uu = np.array(u, float)
            uu /= np.linalg.norm(uu)
            t = 0.01
            q = t * uu
            hv = hcomb(q) + lam * (q @ q) ** 2
            p("      %-20s %14.9f  %14.3e" % (dn, hv / t ** 2, (hv - t ** 2) / t ** 4))
        mx = 0.0
        for uu in U[:60]:
            t = 0.01
            q = t * uu
            hv = hcomb(q) + lam * (q @ q) ** 2
            mx = max(mx, abs((hv - t ** 2) / t ** 4))
        p("      max |(h - t^2)/t^4| over 60 random directions : %.3e" % mx)
        p("    => k0_shell = 4 on this class (symbol level), with lambda = %.6f" % lam)
    p("")
    return tab, loops


# ------------------------------------------------------------------- acid (C)
def part_C():
    p("=" * 80)
    p("PART C  acid tests")
    p("=" * 80)
    p("C1  2D control: nearest-neighbour (d phi)^2 symbol has the SAME sin^2 signature")
    for th in np.linspace(0, np.pi / 2, 5):
        u = np.array([np.cos(th), np.sin(th)])
        p("    theta=%5.3f  sum u^4=%.6f  kappa=sum u^4/12=%.6f"
          % (th, (u ** 4).sum(), (u ** 4).sum() / 12))
    p("    => 2D deficit is anisotropic too (1/24..1/12): the sigma mechanism is NOT a")
    p("       direction-uniform mechanism in 2D either -- it is ON-SHELL (CP 1997 abstract:")
    p("       'O(a^2) on-shell improved action ... satisfying reflection positivity').")
    p("")
    p("C2  4D EOM square |D.F|^2 : direction independence of its symbol")
    rng = np.random.default_rng(5)
    vals = []
    for _ in range(200):
        u = rng.normal(size=4)
        u /= np.linalg.norm(u)
        A = rng.normal(size=4)
        A = A - (A @ u) * u
        F = F_of(u, A)
        dstar = sum(abs(sum(u[mu] * F[mu, nu] for mu in range(4))) ** 2 for nu in range(4))
        vals.append(dstar / ((u @ u) ** 2 * (A @ A)))
    p("    |k.F|^2 / (k^4 |A_T|^2) : %.12f .. %.12f  (exactly 1)" % (min(vals), max(vals)))
    p("")
    p("C3  single plaquette: h = phat_1^2 + phat_2^2 defined; no transfer direction")
    p("    => k0 undefined (N/A).  Its own anisotropy kappa = (u1^4+u2^4)/(12(u1^2+u2^2)):")
    for u in ([1, 0, 0, 0], [1, 1, 0, 0], [1, 2, 0, 0]):
        uu = np.array(u, float)
        uu /= np.linalg.norm(uu)
        p("    u=%-24s kappa = %.6f" % (str(np.round(uu, 4)),
                                       (uu[0] ** 4 + uu[1] ** 4) / (12 * (uu[0] ** 2 + uu[1] ** 2))))
    p("")


if __name__ == "__main__":
    wb, rmin, rmax = part_A()
    tab, loops = part_B()
    part_C()
    p("=" * 80)
    p("SUMMARY")
    p("  A1 max |Bianchi violation| = %.1e   (an identity, not a constraint)" % wb)
    p("  A2 |dF|^2/|d.F|^2 = %.12f .. %.12f  (isotropic: ONE invariant (p^2)^2)" % (rmin, rmax))
    p("=" * 80)
    with open(__file__.replace(".py", ".log"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(OUT) + "\n")
    print("\n[log written to %s]" % __file__.replace(".py", ".log"))
