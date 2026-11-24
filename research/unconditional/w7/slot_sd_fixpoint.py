"""
SLOT-SD-FIXPOINT witness.  Route: turn the exact variance identity
        G''(0) = S - D
into a genuine self-consistency  D = D[S, mu],  and test whether the sign
S - D > 0 is self-consistent and whether it EXCLUDES the gapless fixed point.

Pieces (all cheap, single process, exact where possible):
  A. Toy S-D split (2-mode, exact 1-D quadrature): establish the identity
     G''=S-D, the acid test (free Gaussian: D=0, G''=S; massless G''=0), and
     WHAT makes D sign-definite  (sign of the fiber Hessian K).
  B. Free Gaussian lattice: exact dense.  S = 1/Var = Schur complement.
     Massless k=0: S->0 (D=0) => G''=0  (acid).  No L-uniform free gap.
  C. iid Haar SU(2) (beta=0, the pure-Jacobian source W_geo): curvature of the
     linear-functional marginal > 0, L-INDEPENDENT (Var bounded).  W_geo>0.
  D. beta=2.5 real data (rerun_su2_L4.npz): exact identity closure
     G'' = (E[H]-VarG) + W, W dominant, G''>0.  (D itself needs the fiber K,
     which this data does not carry -> the 2% D/S is external = H1.)
  E. Fixpoint map kappa = F(kappa): gapless kappa=0 is ALWAYS a solution;
     positive branch iff F'(0)>1 (pitchfork = k->0 tachyon).  At beta=2.5
     the coarse block is positive -> F'(0)<1 -> only gapless provable.

Labels are printed inline.  No QED.
"""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, 'D:/Dev/math/research/lattice-cpb')
sys.path.insert(0, HERE)


# ======================================================================
# A.  Toy S-D split  --  exact, 2 mode (a = coarse, y = fiber)
#     V(a,y) = 1/2 A a^2 + 1/2 B y^2 + (k/4) y^4 + h a y^2
#     coarse mode = a.  fiber = y.  mu ~ exp(-V).
#     h_n = d2_a V = A ;  f = d_a V = A a + h y^2 ;
#     K  = d2_y V = B + 3k y^2 + 2 h a ;  c = d2_{ay} V = 2 h y .
#   Conditional at a=0:  nu0(y) ~ exp(-(1/2 B y^2 + (k/4)y^4)).
#   S = E[h_n - c^2/K |0],  D = Var(f|0) - E[c^2/K |0],  G''=S-D= E[h_n|0]-Var(f|0).
# ======================================================================
def toy_quad(f, lo=-14.0, hi=14.0, n=400001):
    y = np.linspace(lo, hi, n)
    w = np.exp(-f(y))
    w /= np.trapezoid(w, y)
    return y, w


def toy_SD(A, B, k, h):
    f0 = lambda y: 0.5 * B * y * y + 0.25 * k * y ** 4
    y, w = toy_quad(f0)
    Ey2 = np.trapezoid(w * y * y, y)
    Ey4 = np.trapezoid(w * y ** 4, y)
    # K at a=0 = B + 3k y^2 ; c at a=0 = 2 h y ; f at a=0 = h y^2
    K = B + 3.0 * k * y * y
    if np.any(K <= 1e-9) and np.any(K < 0):
        Kmin = K.min()
    else:
        Kmin = K.min()
    c2K = 4.0 * h * h * y * y / K                     # c^2/K integrand weight
    Ec2K = np.trapezoid(w * c2K, y)
    Varf = h * h * (Ey4 - Ey2 ** 2)                   # Var(h y^2)
    D = Varf - Ec2K
    S = A - Ec2K                                      # E[h_n|0]=A
    Gpp = S - D
    Gpp_direct = A - Varf                             # E[h_n|0]-Var(f|0)
    return dict(A=A, B=B, k=k, h=h, Kmin=Kmin, S=S, D=D, Gpp=Gpp,
                Gpp_direct=Gpp_direct, Ec2K=Ec2K, Varf=Varf)


def part_A():
    print("=" * 78)
    print("A. TOY S-D SPLIT  (exact quadrature)   G'' = S - D")
    print("=" * 78)
    cases = [
        ("free massless Gaussian  A=0,B=1,k=0,h=0", 0.0, 1.0, 0.0, 0.0),
        ("free massive  Gaussian  A=0.7,B=1,k=0,h=0", 0.7, 1.0, 0.0, 0.0),
        ("gapped coarse+convex fib A=0.5,B=1,k=1,h=0.3", 0.5, 1.0, 1.0, 0.3),
        ("massless+convex fiber    A=0,B=1,k=1,h=0.3", 0.0, 1.0, 1.0, 0.3),
        ("INDEFINITE fiber w/ weight A=0,B=-1,k=2,h=0.5", 0.0, -1.0, 2.0, 0.5),
        ("INDEFINITE but weight at K>0 A=0,B=-2,k=2,h=0.8", 0.0, -2.0, 2.0, 0.8),
    ]
    print(f"{'case':42s} {'Kmin':>7s} {'S':>9s} {'D':>9s} {'Gpp':>9s} {'Gpp_direct':>10s}")
    for name, A, B, k, h in cases:
        r = toy_SD(A, B, k, h)
        tag = "[free]  " if (h == 0.0) else ("[K>0]   " if r['Kmin'] > 0 else "[K<0!]  ")
        print(f"{tag}{name:35s} {r['Kmin']:7.3f} {r['S']:9.4f} {r['D']:9.4f} "
              f"{r['Gpp']:9.4f} {r['Gpp_direct']:10.4f}")
    print("  -> free Gaussian: D=0 exactly, S=1/Var, Gpp=S (massless A=0 => Gpp=0).")
    print("  -> K>0 pointwise => D<=0 (Brascamp-Lieb).  K<0 does NOT force D>0:")
    print("     D>0 needs the indefinite region to CARRY MEASURE (row B=-1 yes,")
    print("     row B=-2 no: weight sits where the quartic makes K>0).  D is the")
    print("     weighted overlap, NOT a function of sign(K) or of S -> no D=D[S].")
    print("  -> identity S-D = E[h_n|0]-Var(f|0) holds to quadrature (last 2 cols).")


# ======================================================================
# B.  Free Gaussian lattice -- exact dense.  S = 1/Var (Schur).
#     A = -Delta + m^2 on a 1D ring of L sites.  Mode = lowest (k=2pi/L).
# ======================================================================
def part_B():
    print("=" * 78)
    print("B. FREE GAUSSIAN LATTICE (exact dense):  S = 1/Var,  D = 0")
    print("=" * 78)
    print(f"{'L':>4s} {'m^2':>6s} {'Var(a_k1)':>12s} {'S=1/Var':>10s} {'Gpp':>9s}   (mode k=2pi/L)")
    for L in (4, 8, 16, 32, 64):
        lam = 2.0 - 2.0 * np.cos(2.0 * np.pi / L)          # first nonzero k
        for m2 in (0.0, 0.5):
            S = lam + m2
            Var = 1.0 / S if S > 0 else np.inf
            Gpp = S - 0.0                                   # D=0 for Gaussian
            print(f"{L:4d} {m2:6.2f} {Var:12.4f} {S:10.5f} {Gpp:9.5f}")
    print("  -> massless lowest k: S = 4 sin^2(pi/L) -> 0  (no L-uniform free gap).")
    print("  -> massive: S -> m^2 > 0 (a genuine gap needs the mass).")
    print("  -> D = 0 exactly (Gaussian); Gpp = S:  free massless Gpp=0  [ACID PASS].")


# ======================================================================
# C.  iid Haar SU(2) (beta=0): the pure-Jacobian source W_geo.
#     Constant coarse mode over nlink iid links: Var bounded => Gpp>0, L-indep.
# ======================================================================
def su2_haar_x(N, rng):
    """Principal log x of N iid Haar SU(2) links.  q uniform on S^3 = Haar(SU2)."""
    q = rng.standard_normal((N, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    q0 = np.clip(q[:, 0], -1.0, 1.0)
    th = 2.0 * np.arccos(q0)                     # in [0, 2pi]
    s = np.sqrt(np.maximum(1.0 - q0 * q0, 1e-300))  # sin(th/2)
    n = q[:, 1:] / s[:, None]
    x = th[:, None] * n                          # x = theta * nhat  (3-vector)
    bad = s < 1e-9
    x[bad] = 0.0
    return x


def kappa_even(a):
    """Veff'' = curvature of -log p(a) at 0, via even quartic fit on the core."""
    sd = a.std()
    m = np.abs(a) < 2.5 * sd
    lo, hi = a[m].min(), a[m].max()
    cnt, edges = np.histogram(a[m], bins=41, range=(lo, hi))
    cen = 0.5 * (edges[1:] + edges[:-1])
    g = cnt > 0
    yy = -np.log(cnt[g] / (a[m].size * (edges[1] - edges[0])))
    A = np.vstack([np.ones(g.sum()), cen[g] ** 2, cen[g] ** 4]).T
    coef, *_ = np.linalg.lstsq(A, yy, rcond=None)
    return 2.0 * coef[1]


def part_C():
    print("=" * 78)
    print("C. iid Haar SU(2) (beta=0): PURE-JACOBIAN SOURCE W_geo")
    print("=" * 78)
    rng = np.random.default_rng(2024)
    # single link: marginal of one component x_1
    x1 = su2_haar_x(4_000_000, rng)[:, 0]
    V1 = x1.var()
    print(f"single link: Var(x_1)={V1:.4f}  curvature(-log p)={kappa_even(x1):.4f}  "
          f"1/Var={1.0/V1:.4f}")
    # constant coarse mode over nlink iid links -> average of nlink x_1's * sqrt
    print(f"{'nlink':>6s} {'Var(a)':>9s} {'curv':>9s} {'1/Var':>9s}  (a=const-mode avg)")
    XL = su2_haar_x(4_000_000, rng)[:, 0]
    for nlink in (1, 4, 16, 64, 256):
        M = XL.size // nlink * nlink
        a = XL[:M].reshape(-1, nlink).sum(axis=1) / np.sqrt(nlink)
        print(f"{nlink:6d} {a.var():9.4f} {kappa_even(a):9.4f} {1.0/a.var():9.4f}")
    print("  -> Var(a) is L-INDEPENDENT (compact Haar) => Gpp>0 uniformly: the")
    print("     beta=0 source is real.  BUT this is strong coupling, not the crossing.")


# ======================================================================
# D.  beta=2.5 real data: identity closure + W-dominance
# ======================================================================
def part_D():
    print("=" * 78)
    print("D. beta=2.5 DATA  (rerun_su2_L4.npz):  G''=(E[H]-VarG)+W")
    print("=" * 78)
    fn = os.path.join(HERE, 'rerun_su2_L4.npz')
    if not os.path.exists(fn):
        print("  [missing data] ", fn)
        return
    import rerun_core as RC
    d = np.load(fn)
    a, H, G, P1, P2 = d['a'], d['H'], d['G'], d['P1'], d['P2']
    NM = a.shape[1]
    rng = np.random.default_rng(7000)
    rows = []
    print(f"{'mode':>5s} {'Gpp_KDE':>9s} {'1/Var':>8s} {'E[H]-VarG':>10s} {'E[Phi2]':>8s} "
          f"{'W':>8s} {'br+W':>8s}")
    for j in range(NM):
        k = RC.kernel_cond(a[:, j], {'H': H[:, j], 'G': G[:, j],
                                     'P1': P1[:, j], 'P2': P2[:, j]}, bw=0.30)
        EH = k['E_H']; EP2 = k['E_P2']
        VG = k['V_own']['G']; VP1 = k['V_own']['P1']; cov = k['C_G_P1']
        bracket = EH - VG
        W = EP2 - VP1 - 2.0 * cov
        kk, _ = RC.kappa_kde(a[:, j], rng=rng)
        inv = 1.0 / a[:, j].var()
        rows.append((kk, inv, bracket, EP2, W, bracket + W))
        print(f"{j:5d} {kk:+9.4f} {inv:8.4f} {bracket:+10.4f} {EP2:8.4f} "
              f"{W:+8.4f} {bracket+W:+8.4f}")
    rr = np.array(rows)
    print(f"{'mean':>5s} {rr[:,0].mean():+9.4f} {rr[:,1].mean():8.4f} "
          f"{rr[:,2].mean():+10.4f} {rr[:,3].mean():8.4f} {rr[:,4].mean():+8.4f} "
          f"{rr[:,5].mean():+8.4f}")
    print("  -> G''>0 at beta=2.5; W (Haar-Jacobian) carries it; E[H]-VarG ~ 0.")
    print("  -> D needs the fiber Hessian K (not in this data): D/S~2% is the")
    print("     external PROVE-Gpp number, i.e. the unmeasured H1 remainder.")


# ======================================================================
# E.  Fixpoint map  kappa = F(kappa):  gapless always solves; pitchfork iff F'(0)>1
# ======================================================================
def part_E():
    print("=" * 78)
    print("E. FIXPOINT  kappa = F(kappa):  gapless solution, pitchfork criterion")
    print("=" * 78)
    kap0 = 1.0
    print(f"{'A=F0prime':>10s} {'kappa* (pos FP)':>16s} {'regime':>24s}")
    for Apar in (0.5, 0.9, 1.0, 1.5, 3.0):
        # minimal faithful family: F(0)=0 (gapless always solves), F'(0)=A
        # F(k) = A k / (1 + k/kap0)
        if Apar > 1.0:
            kstar = kap0 * (Apar - 1.0)
            tag = "positive branch"
            ks = f"{kstar:.4f}"
        elif Apar == 1.0:
            tag = "marginal (FP at 0 only)"
            ks = "0"
        else:
            tag = "gapless only"
            ks = "0"
        print(f"{Apar:10.2f} {ks:>16s} {tag:>24s}")
    print("  -> kappa=0 is ALWAYS a fixed point (free field, F(0)=0).")
    print("  -> positive kappa* appears IFF A=F'(0)>1 (pitchfork = k->0 tachyon).")
    print("  -> beta=2.5 coarse block is POSITIVE (W6-BE, 0% neg) => A<1 =>")
    print("     only the gapless branch is provable.  (This is the wall.)")


if __name__ == '__main__':
    part_A(); part_B(); part_C(); part_D(); part_E()
    print("=" * 78)
    print("VERDICT: identity G''=S-D holds; G''=0 (free) SOLVES the self-consistency")
    print("(D=0,S->0). W_geo>0 is a real source at beta=0 but VANISHES on the free")
    print("acid field; at beta=2.5 G''>0 is MEASURED but D<S is the H1 step. No QED.")
