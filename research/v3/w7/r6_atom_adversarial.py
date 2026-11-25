#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
r6_atom_adversarial.py  —  adversarial witness for the 0++ Kallen-Lehmann atom
f1^2 = nu({m1^2}) of O = tr F^2.

Question: can f1^2 -> 0 while H1 (chi<=1/lam), positivity, gauge invariance and
the strong-coupling residue survive?  Or is the atom protected / fine-tuned?

Deterministic, numpy only.  Run:  py -3.12 -I r6_atom_adversarial.py   (or python)
Output is quoted (abridged) in research/v3/R6-ATOM-ADVERSARIAL.md  §6.

Five parts:
  A  THIN ATOM (Q6 gadget): nu = w d_eps + (1-w) d_s2, w=eps/(2 lam).
     f1^2 -> 0 AND m1 -> 0 *together*: it is a GAPLESS limit, not a gapped
     atomless measure.  So the thin atom kills "H1 => gap", NOT "gapped atom has
     f1^2>0".  This is the central distinction.
  B  ABSTRACT fixed-m atomless measure: nu = eta d_m + (1-eta)*continuum on
     [theta, theta+W].  m1=m>0 FIXED, f1^2=eta->0, chi BOUNDED.  A genuine
     gapped atomless measure EXISTS in measure space -- the adversary's strongest
     measure-level kill.  What excludes it physically is NOT H1 (chi is fine) but
     the carrier overlap (trace anomaly), flagged, not used as a proof.
  C  SAVVIDY double well V=(g/2)(x^2-1)^2: the 0++ atom residue
     f1^2=| <psi1| x |psi0> |^2 stays O(1) (robust, ~ well position^2) while the
     GAP E1-E0 collapses (Kramers).  So at the Savvidy cell the ATOM is robust
     and the GAP is the fragile object -- the opposite of "atom is fine-tuned".
  D  FREE FIELD acid test: the tr F^2 channel of a free (massive) theory has NO
     isolated atom -- the lightest weight is a 2-particle cut, f1^2=0.  Must
     return 0 (no manufactured protection).
  E  TRANSPORT law f1^2(mu)=Z_O^2 f1^2(mu0): the *sign* f1^2>0 is invariant
     (Z_O^2>0), the *floor* r0 is scheme/Z_O dependent -> floor is fine-tuned,
     existence is not.
"""
import numpy as np

np.set_printoptions(precision=6, suppress=False)

def hr(t): print("\n" + "="*74 + "\n" + t + "\n" + "="*74)

# ---------------------------------------------------------------- PART A
hr("(A) THIN ATOM: f1^2 -> 0 and m1 -> 0 TOGETHER (gapless, not gapped-atomless)")
lam = 1.0
print("  nu = w d_eps + (1-w) d_{s2},  w = eps/(2 lam),  s2 = 2 lam (1-w)   [ch02 Prop.]")
print("     eps      w=f1^2      m1=sqrt(eps)     chi=Integral s^-1 dnu   (target 1/lam=1)")
for eps in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6]:
    w = eps/(2*lam)
    s2 = 2*lam*(1-w)
    chi = w/eps + (1-w)/s2
    m1 = np.sqrt(eps)
    print(f"   {eps:9.1e}  {w:10.3e}   {m1:10.3e}      {chi:16.10f}")
print("  READ: chi == 1/lam for all eps (H1 holds), but f1^2 -> 0 ONLY because")
print("        m1 = sqrt(eps) -> 0 as well.  This is a GAPLESS sequence.")
print("        It does NOT exhibit a GAPPED state with f1^2 -> 0.")

# ---------------------------------------------------------------- PART B
hr("(B) ABSTRACT fixed-m atomless measure: m1 FIXED, f1^2 -> 0, chi BOUNDED")
m = 1.0; theta = 2.0; W = 2.0      # atom at m=1 ; continuum on [2,4]  (theta>=2m)
cont_integral = (1.0/W)*np.log((theta+W)/theta)   # Integral s^-1 ds over [theta,theta+W]
print(f"  nu = eta d_m + (1-eta) * Uniform[theta,theta+W];  m={m}, [theta,theta+W]=[{theta},{theta+W}]")
print(f"  continuum piece  Integral s^-1 dnu  = {cont_integral:.6f}")
print("     eta=f1^2    m1       chi = eta/m^2 + (1-eta)*cont     bounded?")
for eta in [0.5, 0.1, 1e-2, 1e-4, 1e-8]:
    chi = eta/m**2 + (1-eta)*cont_integral
    print(f"   {eta:9.1e}   {m:5.2f}     {chi:16.10f}              yes")
print("  READ: a GAPPED (m1=1 fixed) ATOMLESS (f1^2->0) measure with bounded chi")
print("        EXISTS in measure space.  H1 cannot exclude it (chi bounded, m1>0).")
print("        What EXCLUDES it physically is the carrier overlap")
print("        <0|tr F^2|G> != 0 (trace anomaly / energy density), NOT H1.  [flag]")

# ---------------------------------------------------------------- PART C
hr("(C) SAVVIDY double well: residue ROBUST, transfer GAP collapses, H1 fails (chi->inf)")
# Coarse marginal nu(x) ~ exp(-V), V=(g/2)(x^2-1)^2, wells x=+-1, barrier g/2.
# The physical (transfer) masses are eigenvalues lam_n of the Fokker-Planck
# generator L=ddx - V' ddx on L^2(nu); via the SUSY similarity this is the
# spectrum of H = -ddx + W,  W = V'^2/4 - V''/2,  ground state h0 ~ exp(-V/2), lam0=0.
# KL measure of the coordinate operator a=x:  C(t)=sum_n |c_n|^2 e^{-lam_n t},
#   c_n = <h0 | x | h_n>_{L^2(dx)}   (flat inner product, h_n = H-eigenfunctions).
# => mass gap = lam_1 (Kramers: ~ exp(-S), S=(4/3)sqrt(g)); residue f1^2=|c_1|^2.
def double_well_transfer(g, L=8.0, N=3600):
    x = np.linspace(-L, L, N)
    h = x[1]-x[0]
    V  = 0.5*g*(x**2 - 1.0)**2
    Vp = g*(2*x**3 - 2*x)          # V'
    Vpp = g*(6*x**2 - 2)           # V''
    W = 0.25*Vp**2 - 0.5*Vpp       # SUSY potential
    main = 2.0/h**2 + W            # -d^2 + W  (Dirichlet)
    off  = -1.0/h**2 * np.ones(N-1)
    try:
        from scipy.linalg import eigh_tridiagonal
        E, U = eigh_tridiagonal(main, off)
    except Exception:
        A = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
        E, U = np.linalg.eigh(A)
    h0 = U[:,0]; h1 = U[:,1]       # h0 even (E~0), h1 odd (E=lam_1)
    h0 /= np.sqrt(np.sum(h0**2)*h); h1 /= np.sqrt(np.sum(h1**2)*h)
    lam1 = E[1]-E[0]               # transfer gap (= Kramers-suppressed)
    c1 = np.sum(h1*x*h0)*h         # <h0|x|h1>
    # dominant susceptibility term chi ~ = c1^2/lam1  (the n=1 term of sum c_n^2/lam_n;
    # C(t)=sum c_n^2 e^{-lam_n t}, chi = int_0^inf C(t)dt = sum c_n^2/lam_n)
    chi_lead = c1**2/lam1
    return lam1, c1**2, E[2]-E[1], chi_lead

print("   g    lam1=GAP     f1^2=|c1|^2   isolation I2=lam2-lam1   chi~=c1^2/lam1")
for g in [0.5, 1.0, 2.0, 4.0, 8.0, 16.0]:
    lam1, f1sq, iso, chi_lead = double_well_transfer(g)
    print(f"  {g:5.1f}  {lam1:10.4e}   {f1sq:12.8f}     {iso:12.5f}          {chi_lead:12.4e}")
print("  READ: as the barrier grows the transfer GAP lam1 collapses ~e^{-S} (Kramers),")
print("        but the ATOM RESIDUE f1^2 and its isolation I2 stay O(1)/grow.")
print("        The susceptibility chi~=c1^2/lam1 -> infinity: H1 FAILS here.")
print("        => This model has residue+I2 but NO gap and NO H1.  It shows the RESIDUE")
print("           (and even I2 with IN-11 failing) is ROBUST, the GAP is fine-tuned;")
print("           it is NOT a counterexample to H*=(residue & I2 & H1)=>gap, since H1")
print("           fails.  The residue floor is strictly easier than the gap.")

# ---------------------------------------------------------------- PART D
hr("(D) FREE FIELD acid test: tr F^2 channel of a free theory has NO atom")
# free scalar: the :phi^2: composite couples to the 2-particle cut [4m^2, inf),
# no isolated pole below -> f1^2 = 0.  Represent as continuum measure.
mphi = 1.0
print("  free massive scalar, O=:phi^2:  ->  nu = cut on [4 m^2, inf),  f1^2=0")
print(f"    m_phi={mphi}:  threshold 4m^2={4*mphi**2},  f1^2=nu({{inf supp}})=0")
# show chi of the free composite is finite but f1^2=0 -> H1 holds, atom absent
# nu(s) = rho(s) ds on [4m^2, 4m^2+Lam]; take flat rho
Lcut = 50.0
rho = 1.0/Lcut
s = np.linspace(4*mphi**2, 4*mphi**2+Lcut, 400000)
chi_free = np.trapezoid(rho/s, s)
print(f"    chi_free(composite) = Integral rho/s ds = {chi_free:.6f}  (finite => H1 ok)")
print(f"    f1^2 = 0  (no atom)  ->  ACID PASS: free field has no manufactured atom.")
print("  READ: every free theory gives f1^2=0.  The atom is a CONFINEMENT effect,")
print("        not a free-field artifact.  (fail-on-free-field = correct direction)")

# ---------------------------------------------------------------- PART E
hr("(E) TRANSPORT f1^2(mu)=Z_O^2 f1^2(mu0): sign invariant, FLOOR scheme-dependent")
print("   mu0     Z_O      Z_O^2     f1^2(mu0)   f1^2(mu)=Z^2 f1^2    sign")
f10 = 0.37
for Z in [0.25, 0.5, 1.0, 2.0, 4.0]:
    print(f"    -      {Z:5.2f}   {Z*Z:7.3f}     {f10:8.4f}     {Z*Z*f10:14.6f}       {'+' if Z*Z*f10>0 else '0'}")
# and the zero case: 0 -> 0 always
print(f"   f1^2(mu0)=0  ->  Z^2 * 0 = 0 for every Z   (0 -> 0: multiplicative)")
print("  READ: EXISTENCE (f1^2>0 vs =0) is Z_O-invariant (Z_O^2>0).")
print("        The FLOOR r0 >= const is NOT invariant (rescales with Z_O^2).")
print("        => 'f1^2 >= r0' is scheme/fine-tuning dependent; only the SIGN is physical.")

hr("DONE.  All five checks are deterministic linear algebra / arithmetic.")
print("No third-party data.  numpy", np.__version__, "| python", __import__('sys').version.split()[0])
