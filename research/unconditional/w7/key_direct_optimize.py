"""
key_direct_optimize.py  --  witness for research/unconditional/KEY-DIRECT-OPTIMIZE.md
Single small job.  Exact linear algebra / 1D-2D quadrature / tiny toys.  NO lattice MC (heavy).
Run:  py -3.12 key_direct_optimize.py

OPTIMIZATION of the Savvidy signed-minorant construct (KEY-DIRECT-DERIVATION attack I):
  (A) build Phi_total(B) = [signed Savvidy part sqrt(2 g B)] + [W_geo >= 0], attempt Phi_total <= m1
      and decide: does it lower-bound m1 (a floor) or only give G''>0 (H1 = inverse susceptibility)?
  (B) the J^{PC} / channel of the Savvidy negative mode: 1-form magnetic (J=1) or 0++ projection?
  (C) H1 -> H* upgrade via an exact identity (Fredholm / completely-monotone): does it FORCE the atom?

Labels [T]/[C]/[O]/[NG]/[V-here] live in the .md; nothing here is a proof.
Nothing is a fake closure: the Clay wall is the expected honest output.
"""
import numpy as np
from scipy import integrate

np.set_printoptions(precision=5, suppress=True)
rng = np.random.default_rng(0)


def banner(s):
    print("\n" + "=" * 80)
    print(s)
    print("=" * 80)


# ===========================================================================
# PART A1  --  W_geo (pure-Haar marginal curvature): NONABELIAN-ONLY positivity.
#   Single link.  SU(2) flat chart X in su(2), |X| = phi in [0,2pi]; Haar density
#   J(phi) = (sin(phi/2)/(phi/2))^2  (Vandermonde; const irrelevant to curvature).
#   Marginal of  a = <v, X> ,  p(a) = int J(phi) d^3X delta(X.v - a).
#   For U(1): flat chart X in [-pi,pi), J=1  ->  p(a) = const  ->  curvature 0.
# ===========================================================================
banner("PART A1  |  W_geo = pure-Haar marginal curvature : SU(2) > 0, U(1) = 0   [nonabelian-only]")

def p_SU2(a, phi_max=2.0 * np.pi):
    """marginal density of one flat component a=|X.v| of SU(2) log, up to a const.
    d^3X = da * rho drho * dpsi ;  phi = sqrt(a^2+rho^2);  rho drho = phi dphi
    => p(a) ∝ ∫_{phi=|a|}^{phi_max} J(phi) phi dphi = ∫ sin^2(phi/2)/phi dphi ."""
    f = lambda phi: np.sin(phi / 2.0) ** 2 / phi
    val, _ = integrate.quad(f, abs(a), phi_max, limit=400)
    return val

def curvature_from_p(pfun, h=1e-3):
    p0 = pfun(0.0); pp = pfun(h); pm = pfun(-h)
    p2 = (pp - 2 * p0 + pm) / h**2          # p''(0)
    p1 = (pp - pm) / (2 * h)                 # p'(0) (=0 by symmetry)
    Vpp = -p2 / p0 + (p1 / p0) ** 2          # V_eff'' = -(log p)'' at 0
    return Vpp

Wgeo_SU2 = curvature_from_p(p_SU2)
print("   SU(2) single link:  p(0)=%.5f   W_geo = V_eff''(0) = % .5f   (>0)" % (p_SU2(0.0), Wgeo_SU2))
# analytic cross-check: near 0, p'(a)=-sign(a) sin^2(a/2)/a ~ -a/4, p''(0)=-1/4, V_eff''=(1/4)/p(0)
print("   analytic cross-check: (1/4)/p(0) = % .5f   (matches quadrature)" % (0.25 / p_SU2(0.0)))
# U(1): p(a)=const on [-pi,pi) -> all derivatives 0
print("   U(1)  single link:  p(a)=const  =>  W_geo = 0   (proved: flat Haar, J=1)")
print("   => the Haar-Jacobian curvature is STRICTLY POSITIVE only for SU(N); it VANISHES for U(1).")
print("      This is a genuine nonabelian-only POSITIVE invariant (Prékopa sign at beta=0). [V-here]")

# ===========================================================================
# PART A2  --  the signed Savvidy part: negative mode NONABELIAN-ONLY (recap of witness A)
#   E^2 = 2|qB|(n+1/2) + 2 qB s_z ,  s_z in {-1,0,+1}.  U(1): photon neutral q=0 -> none.
#   SU(2): off-diagonal gluons q=+/-1 -> s_z=-1,n=0 mode E^2=-|qB|<0.
# ===========================================================================
banner("PART A2  |  signed Savvidy part sqrt(2 g B): nonabelian-only instability")
def vect_levels(qB, szs=(-1, 0, 1), nmax=0):
    out = []
    for n in range(nmax + 1):
        for s in szs:
            out.append(2.0 * abs(qB) * (n + 0.5) + 2.0 * qB * s)
    return np.array(sorted(out))

E2_U1 = vect_levels(0.0)
E2_SU2 = np.sort(np.concatenate([vect_levels(1.0 * q) for q in (+1.0, -1.0, 0.0)]))
print("   U(1): min E^2 = %.3f  negative modes = %d" % (E2_U1.min(), (E2_U1 < -1e-9).sum()))
print("   SU(2)@gB=1: min E^2 = %.3f  negative modes = %d  (E^2=-gB)" %
      (E2_SU2.min(), (E2_SU2 < -1e-9).sum()))
print("   Phi_Sav(B) = sqrt(2 g B):  B=1 %.3f | B=0.01 %.3f | B=0 %.3f  (monotone -> 0)"
      % (np.sqrt(2 * 1 * 1), np.sqrt(2 * 1 * 0.01), np.sqrt(2 * 1 * 0.0)))

# ===========================================================================
# PART A3  --  Phi_total(B) = sqrt(2 g B) + W_geo  : STRICTLY POSITIVE at B=0, but a CURVATURE.
#   Show Phi_total(0)=W_geo>0 (task's literal ask ACHIEVED) ...
#   ... then the DECISIVE point: any object built this way is G''=1/chi, an inverse
#       susceptibility.  The exact squeeze  f1^2/chi <= m1^2 <= M0/chi  makes G'' an
#       UPPER bound on m1^2 (harmonic mean >= min).  Witness: dnu = sqrt(s) ds has
#       chi<inf, G''>0, yet inf supp nu = 0.  So G''>0  does NOT  floor m1.
# ===========================================================================
banner("PART A3  |  Phi_total(0)=W_geo>0 achieved -- but it is a CURVATURE (G''=1/chi), not a floor")

def moments_from_measure(s, w):
    """s: support points (array), w: weights (density*nodes).  Return M0, chi, Gpp, m1."""
    w = np.asarray(w); s = np.asarray(s)
    M0 = np.sum(w)
    chi = np.sum(w / np.maximum(s, 1e-300))
    m1sq = s[np.argmax(w > 1e-12)]           # smallest s with weight = inf supp (approx)
    return M0, chi, (1.0 / chi if chi > 0 else np.inf), m1sq

print("\n   (i) WITNESS  dnu = sqrt(s) ds on [0,Lam]  (PURE continuum, NO atom, inf supp = 0):")
Lam = 4.0
s = np.linspace(1e-6, Lam, 2_000_000)
w = np.sqrt(s) * (s[1] - s[0])
M0, chi, Gpp, m1sq = moments_from_measure(s, w)
print("       Lam=%.1f :  M0=%.4f  chi=%.4f  G''=1/chi=%.4f (>0!)   f1^2=0   m1=0"
      % (Lam, M0, chi, Gpp))
print("       => G''>0 AND f1^2=0 AND m1=0  coexist.  (curvature positivity is BLIND to a gapless edge)")

print("\n   (ii) THIN-ATOM witness  nu = w*delta_eps + (1-w)*delta_S  (gapped-looking, eps->0):")
for eps, ww, S in [(0.02, 0.02, 4.0), (0.02, 0.10, 4.0), (0.001, 0.02, 4.0)]:
    M0 = ww + (1 - ww)
    chi = ww / eps + (1 - ww) / S
    Gpp = 1.0 / chi
    print("       eps=%.4f w=%.2f S=%.1f :  G''=%.4f  chi=%.4f  m1^2=eps=%.4f (->0 as eps->0)"
          % (eps, ww, S, Gpp, chi, eps))

print("\n   (iii) SQUEEZE direction: with M0 normalised,  G'' = 1/chi >= m1^2/M0  (harmonic mean >= min).")
print("         So G'' is an UPPER bound on m1^2 -- the OPPOSITE of a floor.  Phi_total<=m1 is NOT a floor.")
print("         [T] The moment  chi=int s^-1 dnu  is a residue-weighted harmonic mean; it CANNOT")
print("             see whether the edge carries an atom (f1^2>0) or is a.c. (f1^2=0).")

# ===========================================================================
# PART B  --  CHANNEL of the Savvidy negative mode: 1-form magnetic (J=1) vs 0++ (A1g).
#   1. COLOR: eigenvector of Hessian in the ADJOINT -> chromo-charged -> ZERO overlap with any
#      color-singlet operator (tr F^2 is a singlet).  [T]
#   2. SPIN/HELICITY: negative mode is the s_z=-1 TRANSVERSE magnetic mode of a massless
#      spin-1 field -> physical helicity +/-1 -> single-particle J>=1 -> NO J=0.  [T]
#   3. BAND: LLL is degenerate (guiding center); J_z = m + s_z spans all integers, but each
#      member is a single adjoint transverse vector.  J_z=0 member exists; J=0 member does NOT.
#   4. TWO-GLUON 0++: minimal color-singlet J=0 needs TWO gluons -> support = cut on [2 m_eff, inf),
#      NOT an isolated atom.  Threshold 2 m_eff = 2 sqrt(g B) > 0.
# ===========================================================================
banner("PART B  |  CHANNEL: the negative mode is a colored 1-form magnetic (J=1) object, not 0++")
# (B.1) color: ad(T3) eigenvalues on su(2) -> the negative-mode gluons are charged (+/-1)
f = np.zeros((3, 3, 3))
for a in range(3):
    for b in range(3):
        for c in range(3):
            if (a, b, c) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]: f[a, b, c] = 1.0
            if (a, b, c) in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]: f[a, b, c] = -1.0
adT3 = f[2]
print("   ad(T3) eigenvalues:", np.linalg.eigvals(adT3),
      " -> negative-mode gluons carry charge q=+/-1 (ADJOINT, chromo-charged).")
print("   tr F^2 (the 0++ atom operator) is a COLOR SINGLET -> <singlet|adjoint>=0 (Schur). [T]")

# (B.2) spin/helicity: s_z of the unstable mode and single-particle J minimum
sz_neg = -1                       # the unstable level is s_z=-1 (spin anti-parallel to B)
helicities = (+1, -1)             # physical polarisations of a massless spin-1
print("   unstable mode: n=0, s_z=%+d  (TRANSVERSE magnetic, spin anti-parallel to B)." % sz_neg)
print("   physical helicities of a massless spin-1: %s  => one-particle J >= |lambda| = 1  => NO J=0. [T]"
      % (str(helicities),))

# (B.3) LLL band: J_z = m + s_z  (symmetric gauge), m=0,1,2,... ; show it spans all J_z
print("\n   LLL guiding-centre band, J_z = m + s_z  (s_z=-1):")
band = [(m, m + sz_neg) for m in range(0, 5)]
for m, Jz in band:
    Jmin = abs(sz_neg)            # minimal total J consistent with helicity/s_z
    print("      m=%d :  J_z=%+d  (member of the band;  still J>=%d, NOT J=0)" % (m, Jz, Jmin))
print("   => the band contains a J_z=0 member (m=1) but that state has J>=1, NOT J=0 (J_z=0 != J=0).")

# (B.4) two-gluon 0++ threshold: color-singlet J=0 = TWO gluons -> continuum from 2 m_eff
print("\n   two-gluon 0++ continuum: minimal color-singlet J=0 state = (gluon)+(gluon),")
for gB in (0.5, 1.0, 2.0):
    m_eff = np.sqrt(gB)
    print("      gB=%.2f :  m_eff=sqrt(gB)=%.3f  2*gluon threshold = %.3f   "
          "(0++ support = [%.3f, inf) : a CUT, no atom below)" % (gB, m_eff, 2 * m_eff, 2 * m_eff))
print("   => the negative mode can reach 0++ ONLY as a two-gluon CUT above 2 m_eff (>0).")
print("      An ISOLATED 0++ atom (pole BELOW the 2-gluon cut) is NOT produced by it.")
print("   => and at the PHYSICAL vacuum <B>=0 the operator is -D^2 >= 0 : the negative mode")
print("      DOES NOT EXIST.  The 'channel bridge' evaporates exactly where the gap is measured. [T]")

# ===========================================================================
# PART C  --  H1 -> H* upgrade test: does G''>0 + an exact identity FORCE f1^2>0 ?
#   C.1: two KL measures with the SAME chi (=same G'') but DIFFERENT f1^2 (atom vs a.c.).
#        => any functional depending only on chi/G'' cannot determine the atom. [T]
#   C.2: Fredholm toy  D_O(z)=det(1 - z T_O): an atom eigenvalue => ISOLATED zero of D_O;
#        a pure cut => NO zero.  D_O READS the atom -- but locating the zero = the spectrum = answer.
# ===========================================================================
banner("PART C  |  H1 -> H* : exact identities READ the atom, they do not FORCE it")

# (C.1) same chi, different f1^2
print("   (C.1) two measures, matched chi (hence matched G''=1/chi), DIFFERENT edge atom f1^2:")
chi_target = 0.75
# (a) atom + continuum: nu = f delta_{m1} + rest; tune rest so chi=chi_target
m1sq, S = 1.0, 9.0
# chi = f/m1sq + (1-f)/S = chi_target  (per unit M0, with M0=1 continuum part density at S)
# solve f/m1 + (1-f)/S = chi_target  -> f(1/m1 - 1/S) = chi_target - 1/S
f_atom = (chi_target - 1.0 / S) / (1.0 / m1sq - 1.0 / S)
# (b) pure a.c. sqrt(s) on [0,Lam]: chi = 2 sqrt(Lam) = chi_target -> Lam = (chi_target/2)^2
Lam2 = (chi_target / 2.0) ** 2
print("       (a) atom+cont:  f1^2 = %.4f (>0)   m1^2=%.2f   chi=%.4f" %
      (f_atom, m1sq, f_atom / m1sq + (1 - f_atom) / S))
print("       (b) sqrt(s) ds on [0,%.3f]: f1^2 = 0   m1^2=0   chi=%.4f" % (Lam2, 2 * np.sqrt(Lam2)))
print("       => SAME chi / SAME G'' ,  f1^2 = %.4f vs 0 ,  m1^2 = %.2f vs 0 ." % (f_atom, m1sq))
print("          A functional of G'' alone CANNOT output the atom or the gap.  (H1-not-H*, [T])")

# (C.2) Fredholm toy: transfer operator spectrum, atom vs cut
print("\n   (C.2) Fredholm determinant toy over a discrete transfer operator T_O (spectrum in (0,1)):")
# toy cut: T has continuous spectrum on [c,1);  toy atom: a single eigenvalue mu below the cut.
def D_cut(z, c=0.5, n=4000):
    lam = np.linspace(c, 1.0, n)             # continuous cut
    return np.prod(1.0 - z * lam)
def D_atom(z, mu=0.2, c=0.5, n=4000):
    lam = np.linspace(c, 1.0, n)
    return (1.0 - z * mu) * np.prod(1.0 - z * lam)
zs = np.linspace(0.0, 0.999, 4000)
Dc = np.array([D_cut(z) for z in zs])
Datom = np.array([D_atom(z) for z in zs])
# zero of D corresponds to eigenvalue 1/z of T.  Atom mu=0.2 -> zero at z=1/mu=5 (OUTSIDE our range 0..1)
# show: inside z in (0,1) the atom contributes an extra factor (1 - z*mu) that vanishes at z=1/mu.
print("       cut-only   : D_cut(z) has NO zero for z in (0,1) (eigenvalues only on the cut [0.5,1]).")
print("       atom+cut   : D_atom(z) has an isolated zero at z = 1/mu = %.2f (mu=0.2)." % (1 / 0.2))
print("       => D_O READS the atom (an isolated zero appears).  BUT its location z1=1/mu = e^{+a m1}")
print("          IS the mass m1 -- locating it uniformly in a is exactly the transport / H1 step.")
print("          The identity does NOT force the atom; it relocates the unknown.  (READ != FORCE, [NG])")

banner("SUMMARY (no QED)")
print(" A: Phi_total(B)=sqrt(2gB)+W_geo is STRICTLY POSITIVE at B=0 (W_geo>0 for SU(N), =0 for U(1)).")
print("    It certifies G'' >= W_geo > 0  (H1).  But G''=1/chi is an INVERSE SUSCEPTIBILITY:")
print("    (i) squeeze puts G'' ABOVE m1^2 (upper bound); (ii) witness dnu=sqrt(s)ds has G''>0,m1=0;")
print("    (iii) thin-atom has G''>0, m1->0.  So it is H1, NOT a floor.  [T]/[V-here]")
print(" B: negative mode = colored (adjoint) transverse-magnetic, s_z=-1, J>=1 single gluon.")
print("    ZERO 0++ projection at one-particle level (color + helicity).  0++ needs TWO gluons ->")
print("    a CUT on [2 m_eff,inf), not an atom; and it does not exist at <B>=0.  WRONG CHANNEL. [T]")
print(" C: H1 + exact identity does NOT give f1^2>0: same chi, different atom; Fredholm READS not FORCES.")
print(" => VERDICT: no optimization yields a 0++ atom floor. Structural death = WRONG-CHANNEL (adjoint")
print("    1-form, not 0++)  AND  H1-not-H*  AND  width-not-floor  AND  no-small-parameter.  Clay 0%.")
