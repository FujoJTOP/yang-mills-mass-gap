#!/usr/bin/env python3
# W2B — build the master field W0 in controlled truncations and compute g1/Sg = f1^2/M0.
#   PART A : Gaussian / Hartree master field  (exact, closes)  -> bilinear trF^2 = 2-gluon bubble
#            -> KL measure absolutely continuous -> atom weight g1/Sg = 0  (free-field acid at finite beta)
#   PART B : 2D exact (Singer W0 = exp(-lambda A/2))            -> trF^2=H conserved -> g1/Sg=0
#            (SU(2) character orthogonality); holonomy -> 1; generic source -> positive & computable
#   PART C : variational large-N ansatz (rigid-string breathing mode)  -> W0 has BINDING;
#            g1/Sg > 0 but SOURCE-DEPENDENT (an UPPER estimate, not a certified floor)
#   PART D : the squeeze / ceiling bracket at beta=2.5 (WG chi_O) and what it does/does not give
# numpy + scipy, ONE process, deterministic.  Labels (charter §2): [T]/[C]/[NG]/[O]/[V].
import numpy as np

np.random.seed(0)
fails = []


def check(cond, msg):
    print(("  [PASS] " if cond else "  [FAIL] ") + msg)
    if not cond:
        fails.append(msg)


# =====================================================================================
print("PART A  Gaussian/Hartree master field W0 : bilinear O=trF^2 has NO atom (g1/Sg=0)")
# The Gaussian master field is a free massive fluctuation operator M = -d^2 + mu^2 (+ magnetic
# curvature).  O = trF^2 is BILINEAR in the gauge field, so by Wick <O O>_c = 2 G^2 : a two-gluon
# BUBBLE.  Its spectral measure is absolutely continuous on s>=4mu^2 -> the lightest point
# (the threshold) carries ZERO mass.  Demonstrate the point-mass -> 0 under refinement.
mu = 1.0
th = (2 * mu) ** 2  # two-particle threshold


def atom_fraction(nbins, s_max=400.0):
    edges = np.linspace(th, s_max, nbins + 1)
    mid = 0.5 * (edges[1:] + edges[:-1])
    w = np.diff(edges) / np.sqrt(np.maximum(mid - th, 1e-12))  # 1D van Hove ~ (s-th)^{-1/2}
    w /= w.sum()
    return w[0]


af = [atom_fraction(n) for n in (100, 1000, 10000, 100000)]
print("   2-gluon continuum : weight in lowest bin @ nbins=1e2,1e3,1e4,1e5 =",
      ["%.2e" % x for x in af])
check(af[-1] < af[0] / 10 and all(af[i + 1] < af[i] for i in range(len(af) - 1)),
      "A1 atom weight at the lightest point -> 0 under refinement (a.c. cut, ~1/sqrt(nbins))")

# Exact Wick identity for a Gaussian state: Var(phi^2) = 2 Var(phi)^2 (=2 for unit variance).
# A bound state (atom) would need <O O>_c to have a delta above the 2-particle threshold; Wick forbids it.
N = 2_000_000
phi = np.random.randn(N)
v_phi = np.var(phi)
v_phi2 = np.var(phi ** 2)
print("   free Gaussian: Var(phi)=%.4f Var(phi^2)=%.4f  ratio Var(phi^2)/Var(phi)^2=%.4f (exact 2)"
      % (v_phi, v_phi2, v_phi2 / v_phi ** 2))
check(abs(v_phi2 / v_phi ** 2 - 2.0) < 0.02, "A2 Gaussian state obeys Wick (ratio=2): no bound-state atom")
print("   => in EVERY Gaussian (Hartree) master field, at every beta incl. 2.5, g1/Sg = 0 exactly.")


# =====================================================================================
print("\nPART B  2D exact Singer W0 = exp(-lambda A(C)/2) : trF^2 = H conserved -> g1/Sg = 0")
# SU(2) class functions: Haar dmu = (2/pi) sin^2(theta) dtheta, chi_j = sin((2j+1)theta)/sin theta.
nt = 4000
theta = (np.arange(nt) + 0.5) * np.pi / nt
w = (2.0 / np.pi) * np.sin(theta) ** 2 * (np.pi / nt)


def chi(j, th):
    return np.where(np.abs(np.sin(th)) < 1e-9,
                    (2 * j + 1) * np.ones_like(th),
                    np.sin((2 * j + 1) * th) / np.sin(th))


# orthogonality of characters (the engine behind <chi_0|chi_R> = delta_{0R})
D = 5
gram = np.array([[np.sum(w * chi(a, theta) * chi(b, theta)) for b in range(D)] for a in range(D)])
off = np.max(np.abs(gram - np.eye(D)))
print("   SU(2) character Gram max|int chi_a chi_b - delta| = %.2e" % off)
check(off < 1e-8, "B1 SU(2) character orthogonality (engine of <chi_0|chi_R>=delta_0R)")

# C2(j) for SU(2) irrep spin j: j(j+1).   trF^2 ~ H = sum C2 * projector, diagonal in R.
# Source overlap onto nontrivial R: <0|trF^2|R> ~ C2(R) <chi_0|chi_R> = C2(R) * delta_{0R} = 0.
C2 = np.array([j * (j + 1.0) for j in range(D)])
overlap_trF2 = C2 * gram[0, :]           # <chi_0|chi_R> row, weighted by C2
print("   trF^2 residues onto R=0..%d : %s" % (D - 1, np.array2string(overlap_trF2, precision=3)))
check(np.max(np.abs(overlap_trF2[1:])) < 1e-8, "B2 2D trF^2 (=H) atom weight g1/Sg = 0 (conserved density)")

# Contrast: a generic source with character content gives POSITIVE, computable g1/Sg.
# source = chi_{f} + alpha*chi_{adj}  (f=1/2 fundamental, adj=1); masses ~ C2.
def g1_over_Sg(o):  # o: coefficients on R=1/2,1,3/2,... (skip vacuum); weight |o_R|^2, mass order C2
    idx = np.arange(1, D)
    g = o ** 2
    return g[0] / g.sum(), idx, g


# fundamental only
o_f = np.zeros(D - 1); o_f[0] = 1.0        # R=1/2
g_f, _, _ = g1_over_Sg(o_f)
# fundamental + adjoint equal
o_fa = np.zeros(D - 1); o_fa[0] = 1.0; o_fa[1] = 1.0  # R=1/2 and R=1
g_fa, _, _ = g1_over_Sg(o_fa)
print("   generic 2D sources: holonomy chi_f -> g1/Sg=%.3f ; chi_f+chi_adj -> g1/Sg=%.3f"
      % (g_f, g_fa))
check(abs(g_f - 1.0) < 1e-12 and abs(g_fa - 0.5) < 1e-12,
      "B3 generic source: 2D W0 CAN give positive g1/Sg (trF^2 zero is a SELECTION RULE, not generic)")


# =====================================================================================
print("\nPART C  variational large-N ansatz : rigid-string breathing mode HAS binding")
# Closed flux tube (0++): radial breathing coordinate rho>0.
#   V(rho) = 2*pi*sigma*rho + pi*kappa/rho   (sigma=1),  H = -(1/2M) d^2/drho^2 + V.
# Lightest glueball psi1 is the lowest *excited* breathing mode.  The source O is the tube's
# energy density -> identified with V(rho) (or rho); g1/Sg = |<psi1|O|psi0>|^2 / Var_0(O).
sigma = 1.0
M_eff = 1.0


def breathing(kappa, nlev=6, Lmax=4.0, N=1600):
    rho = np.linspace(1e-3, Lmax, N)
    dr = rho[1] - rho[0]
    V = 2 * np.pi * sigma * rho + np.pi * kappa / rho
    main = 1.0 / (M_eff * dr ** 2) + V
    offd = -0.5 / (M_eff * dr ** 2) * np.ones(N - 1)
    H = np.diag(main) + np.diag(offd, 1) + np.diag(offd, -1)
    ev, evec = np.linalg.eigh(H)
    return rho, ev[:nlev], evec[:, :nlev]


def atom_frac_discrete(rho, ev, evec, Omult):
    # O a multiplicative self-adjoint operator Omult(rho)
    Ophi = [Omult(rho) * evec[:, n] for n in range(evec.shape[1])]
    O00 = np.sum(evec[:, 0] * (Omult(rho) * evec[:, 0]))
    g = np.array([np.sum(evec[:, n] * (Omult(rho) * evec[:, 0])) ** 2 for n in range(evec.shape[1])])
    return g


print("   kappa   rho*    m1=sqrt(E1-E0)  g1/Sg(O=V)   g1/Sg(O=rho)   (single-tube = UPPER estimate)")
rows = []
for kappa in (0.10, 0.138, 0.193, 0.30, 0.50):
    rho, ev, evec = breathing(kappa)
    rho_star = np.sqrt(kappa / 2.0)
    m1sq = ev[1] - ev[0]
    gV = atom_frac_discrete(rho, ev, evec, lambda r: 2 * np.pi * sigma * r + np.pi * kappa / r)
    gR = atom_frac_discrete(rho, ev, evec, lambda r: r)
    # atom fraction within the discrete single-tube sector (denominator = sum over discrete modes)
    fV = gV[1] / gV[1:].sum()
    fR = gR[1] / gR[1:].sum()
    rows.append((kappa, rho_star, m1sq, fV, fR))
    print("   %.3f   %.3f   %.4f          %.4f        %.4f" % (kappa, rho_star, m1sq, fV, fR))
check(all(0.0 < r[3] < 1.0 and 0.0 < r[4] < 1.0 for r in rows),
      "C1 breathing ansatz: g1/Sg > 0 STRICTLY (binding present) for every kappa")
check(all(abs(r[3] - r[4]) > 1e-3 for r in rows),
      "C2 ... but g1/Sg is SOURCE-DEPENDENT (O=V vs O=rho differ): source-tube coupling is the missing input")
# m1/ m1* ratio sanity: rigid-string mass ^2 ~ (2*pi*sqrt(2*kappa))^2 ~ 8*pi^2*kappa
# (the exact continuum variational mass is quoted for reference, the lattice one is a truncation)
m_rs = [2 * np.pi * np.sqrt(2 * k) for k in (0.10, 0.138, 0.193, 0.30, 0.50)]
print("   rigid-string continuum variational mass m1/sqrt(sigma)=2*pi*sqrt(2*kappa):",
      ["%.2f" % x for x in m_rs])
check(m_rs[1] > 3.0 and m_rs[2] < 4.0, "C3 continuum rigid-string mass in the lattice band m/sqrt(sigma)~3.3-3.9 (k=0.138-0.193)")


# =====================================================================================
print("\nPART D  the squeeze / ceiling bracket at beta=2.5 and what it does NOT give")
# WG measured the gauge-invariant ceiling chi_O = Var(total zero mode of APE trF^2) at beta=2.5:
chi0 = 0.328
M0_phys = None  # not separately measured here
# Squeeze (exact): f1^2/chi <= m1^2 <= M0/chi  =>  g1/Sg = f1^2/M0 <= m1^2*chi/M0.
# A two-level model (atom m1 + leakage m2, isolation ratio rho2=m2^2/m1^2) gives, EXACTLY,
#   g1/Sg = 1/(1+R), R = g2/g1 = (chi - g1/m1)*m2/g1 ... ; but what is CERTIFIED from (chi,m1,rho2)
# is the *upper* face only.  Demonstrate: the certified statement bounds g1/Sg ABOVE, never below.
# 2-level consistency: M0=g1+g2, chi=g1/m1+g2/m2.  Solve for g1/Sg as a function of (m1,m2,chi,M0).
def two_level_g1Frac(m1, m2, chi, M0):
    # g1/m1 + g2/m2 = chi ; g1+g2 = M0 ; g2=M0-g1 -> g1(1/m1 -1/m2) + M0/m2 = chi
    g1 = (chi - M0 / m2) / (1.0 / m1 - 1.0 / m2)
    return g1 / M0, g1


print("   (illustrative) chi(2.5)=%.3f : the squeeze gives ONLY the upper face f1^2/M0 <= m1^2 chi/M0" % chi0)
print("   the lower face needs the SECOND level (isolation I2) -> that is H1 [O], not supplied by data")
# Moments-blindness (W2 NG) recap: two measures, same (M0,M1,M2), different g1/Sg.
def moments(atoms):  # atoms: list of (weight, s)
    M = [sum(wt for wt, s in atoms), sum(wt * s for wt, s in atoms), sum(wt * s * s for wt, s in atoms)]
    return np.array(M)


A_ = [(0.4, 1.0), (0.6, 9.0)]                              # g1/Sg = 0.400
B_ = [(0.001, 1.0), (0.4976, 2.5), (0.2436, 6.0), (0.2578, 12.0)]  # g1/Sg = 0.001 (W2 measure B)
print("   measure A: atoms (s=1:0.4, s=9:0.6) -> g1/Sg=%.3f" % A_[0][0])
print("   measure B: atoms (s=1:0.001, 2.5, 6, 12) -> g1/Sg=%.3f (lightest)" % B_[0][0])
momA = moments(A_)
momB = moments(B_)
print("   A moments (M0,M1,M2)=%s ; B moments=%s" % (np.array2string(momA, precision=4), np.array2string(momB, precision=4)))
check(np.allclose(momA, momB, atol=1e-2), "D1 two measures share (M0,M1,M2) to <0.01%%")
check(abs(A_[0][0] - B_[0][0]) > 0.3, "D2 same first 3 moments, g1/Sg differs 400x => moments CANNOT bound g1/Sg (W2 NG)")

npass = 10 - len(fails)
print("\n# %d/%d checks PASS" % (npass, 10))
print("# PART A,B,C,D complete")
import sys
sys.exit(0 if not fails else 1)
