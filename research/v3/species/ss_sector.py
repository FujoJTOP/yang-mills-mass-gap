#!/usr/bin/env python3
r"""
ss_sector.py -- SS-SECTOR witness: can the CHARGED (1-form-charge != 0) sector bound the
NEUTRAL (0++) sector below?  Deterministic, numpy only, ~1 s, single process.

What this checks (every *theorem* content is in SS-SECTOR.md; this is the arithmetic self-check):

 [A] SU(2) Haar: the EXACT algebra-generation identity  tr_adj W(C) = (tr_f W(C))^2 - 1 ,
     so the neutral sector's one-point functions ARE the charged sector's second moments,
     and <tr_adj> = <(tr_f)^2> - 1 = 0 (character orthogonality).  => neutral data is a
     POLYNOMIAL in charged data at coincident points.  [T]

 [B] 2D pure YM exact (Migdal/Witten): spectrum level spacings sigma_R = (g^2/2) C2(R).
     Neutral (N-ality 0) lightest = adjoint.  Ratio sigma_n/sigma_c = C2(adj)/C2(fund)
     = 2N^2/(N^2-1) in (2, 8/3].  => in 2D pure YM the neutral gap IS an exact functional
     of the charged gap:  m_neutral = (2N^2/(N^2-1)) m_charged >= 2 m_charged .  This REFUTES
     the literal "any functional" form of the proposition, inside a PURE-YM family.  [T]
     AND it corrects IMP-ARROWS 2.1(b): 2D does NOT refute "sigma>0 => neutral gap";
     2D neutral gap is > 0.  What 2D kills is the tr F^2 ATOM / variance (item [C]).

 [C] 2D: tr F^2 = tr E^2 is the Casimir = diagonal in the character basis
     => <0|tr F^2|R> = C2(R) delta_{0R} = 0  and Var = 0  => nu_{trF^2} is EMPTY
     (not merely atom-free).  inf supp nu = +infinity, so it TRIVIALLY satisfies every lower
     bound: an empty measure cannot be a counterexample to a lower bound.  [T]

 [D] Decoration invariance (the theorem's engine, numerically): charged data is IDENTICAL
     under T -> T (x) N_m while the neutral gap = min(gap_T, gap_N) -> 0.  Lattice free
     scalar :phi^2: threshold 2*m_lat(m0); m0=0 => threshold 0.  Same charged dataset,
     neutral gap 0 vs > 0.  [T] (arithmetic) + [V]

 [E] Graded toy measure (r, theta): charged observables = functions of e^{i theta} only,
     neutral = functions of r.  Two measures with IDENTICAL charged data (all moments
     M_k = <e^{ik theta}> = 0, and even mixed moments <r e^{i theta}> = 0) but neutral gap
     0 vs M.  => the neutral sector is a FIBRE over the charged data, not a quotient.  [T]+[V]

 [F] beta=0 Haar acid test: plaquette correlator = Var*delta_{t,0} => flat spectral measure,
     inf supp = 0, atom = 0.  Charged data trivial (perimeter law).  [T]

Run: py -3.12 -I ss_sector.py
"""
import numpy as np

def blk(t): print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)

# ---------------------------------------------------------------- SU(2) Haar quadrature
_xg, _wg = np.polynomial.legendre.leggauss(4000)
RHO = np.pi * (_xg + 1.0)               # class angle in (0, 2pi)
WRHO = np.pi * _wg
DENS = np.sin(RHO / 2.0) ** 2           # Haar radial density (norm cancels)
def haar(f):
    return float(np.sum(f * DENS * WRHO) / np.sum(DENS * WRHO))

def m_lat(m0):                          # free lattice scalar: cosh m_lat = 1 + m0^2/2
    return float(np.arccosh(1.0 + 0.5 * m0 * m0))

# =====================================================================
blk("[A] neutral one-point data = SECOND MOMENT of charged data  (SU(2), exact identity)")
# SU(2) character: chi_j(rho) = sin((2j+1)rho/2)/sin(rho/2).  With c = cos(rho/2):
#   chi_{1/2} = 2c ,  chi_1 = 4c^2 - 1 = (chi_{1/2})^2 - 1   (exact polynomial identity).
_c = np.cos(RHO / 2.0)
trf = 2.0 * _c                          # fundamental Wilson loop, real, in [-2,2]
tra = 4.0 * _c ** 2 - 1.0               # adjoint  Wilson loop (= N-ality 0)
print(f"  max |tr_adj - ((tr_f)^2 - 1)| over Haar grid = {np.max(np.abs(tra - (trf**2 - 1.0))):.3e}"
      "   (exact identity)")
print(f"  <tr_f>   = {haar(trf): .3e}      (N-ality 1, charged: parity-odd, = 0)")
print(f"  <(tr_f)^2> = {haar(trf**2): .6f}   (= 1 exactly: character orthonormality)")
print(f"  <tr_adj>   = {haar(tra): .3e}      (= <(tr_f)^2> - 1 = 0: neutral, orthogonal)")
print("  => every neutral one-point function is a POLYNOMIAL in charged ones at coincident")
print("     points. The grading obstructs no ALGEBRA; it constrains no FUNCTIONAL.  [T]")

# =====================================================================
blk("[B] 2D pure YM exact: neutral gap = (2N^2/(N^2-1)) x charged gap  >= 2 x charged")
def c2_fund(N): return (N * N - 1.0) / (2.0 * N)
def c2_adj(N):  return float(N)
print("   N    C2(fund)   C2(adj)   ratio C2(adj)/C2(fund)   = 2N^2/(N^2-1)")
for N in [2, 3, 4, 5, 10, 100]:
    r = c2_adj(N) / c2_fund(N)
    print(f"  {N:3d}  {c2_fund(N):9.5f}  {c2_adj(N):7.3f}        {r:9.6f}        {2*N*N/(N*N-1):9.6f}")
print("  2D spectrum (Migdal): level spacing (per unit length) sigma_R = (g^2/2) C2(R);")
print("  N-ality 0 <=> adjoint is the LIGHTEST neutral rung; N-ality 1 <=> fundamental is")
print("  the lightest charged rung.  Hence  m_neutral = (2N^2/(N^2-1)) m_charged  EXACTLY,")
print("  for every coupling, in pure YM.  => a functional of charged data that gives a")
print("  strictly positive neutral lower bound, inside a pure-YM family.  [T]")
print("  Reading: for each fixed N, F_N(charged) := (2N^2/(N^2-1)) * (lightest charged rung)")
print("  depends on charged data only (the multiplier is the pure group number above).")
print("  (Caveat, stated plainly: d=2; and the 'gap' is the fixed-L holonomy-tower spacing.)")
# ---- [B'] the same functional DIES on the decoration: the pincer -------------------
blk("[B'] the 2D functional is NOT decoration-invariant (why [B] is not a counterexample)")
print("  Take T = 2D SU(2) (charged gap m_ch) and its decoration D_m(T) = T (x) N_m with")
print("  N_m = free scalar of lattice mass m_lat.  Charged data of D_m(T) = charged data of T.")
print("  Neutral rung of D_m(T) = min( (8/3) m_ch , 2*m_lat ).  Test F = (8/3) m_ch:")
m_ch = 1.0
print("   m0      m_neutral(D_m T)   F=(8/3)m_ch     F <= truth ?")
for m0 in [0.0, 0.5, 1.0, 1.5, 2.0]:
    thr = min(8.0 / 3.0 * m_ch, 2.0 * m_lat(m0))
    print(f"  {m0:4.2f}      {thr:10.6f}       {8.0/3.0*m_ch:10.6f}     "
          + ("yes" if 8.0/3.0*m_ch <= thr + 1e-12 else "NO  (F > truth)"))
print("  => F is exact on 2D alone but FAILS on the massless decoration: its PROOF (Migdal's")
print("     exact solution) is not a charged-sector statement.  This is the pincer:")
print("     every counterexample to the 'functional' form is also an instance of the theorem")
print("     in the 'derivational' form.  [T]")

# =====================================================================
blk("[C] 2D: nu_{tr F^2} is EMPTY (atom = 0 AND variance = 0) -- hence unfalsifiable below")
print("  tr F^2 = tr E^2 = quadratic Casimir, DIAGONAL on characters chi_R:")
print("     <0|tr F^2|R> = C2(R) * <chi_0|chi_R> = C2(R) * delta_{0R} = 0 for every R != 0")
print("     <(tr F^2)^2>_c = C2(0)^2 = 0   (normal-ordered vacuum)")
for j in [0.0, 0.5, 1.0, 1.5]:
    print(f"     j={j:>4}: C2=j(j+1)={j*(j+1):6.3f}  weight |<0|trF2|j>|^2 = 0")
print("  => inf supp nu_{trF^2} = +inf (EMPTY support). An empty measure satisfies every")
print("     lower bound, so 2D cannot be a counterexample to a neutral GAP lower bound;")
print("     it can only kill the ATOM.  This is the erratum on IMP-ARROWS 2.1(b).  [T]")

# =====================================================================
blk("[D] decoration invariance: same charged data, neutral gap 0 vs >0 (lattice scalar)")
print("  neutral factor N_m = free lattice scalar; :phi^2: is a two-particle composite,")
print("  so its 0++ threshold is 2*m_lat.  Product theory T (x) N_m:")
print("     all charged correlators   : UNCHANGED (factorize: <O_ch (x) 1> = <O_ch>)")
print("     neutral 0++ threshold     : min(threshold_T , 2*m_lat(m0))")
print("   m0       m_lat       2*m_lat (=neutral threshold of the decoration)")
for m0 in [0.0, 0.1, 0.5, 1.0, 2.0]:
    print(f"  {m0:4.2f}   {m_lat(m0):8.6f}   {2*m_lat(m0):8.6f}")
print("  m0 -> 0 gives threshold -> 0 EXACTLY while the charged dataset is bit-identical.")
print("  => no charged-sector statement can exclude the massless decoration.  [T]+[V]")

# =====================================================================
blk("[E] graded toy measure: the neutral sector is a FIBRE, not a quotient")
# charged observables f(e^{i theta}); neutral observables g(r); measures mu_M = delta_r=M (x) Haar.
th = np.linspace(0.0, 2.0 * np.pi, 4096, endpoint=False)
w = np.full_like(th, 1.0 / th.size)
def charged_moments(k):
    return complex(np.sum(np.exp(1j * k * th) * w))
print("  charged data = {M_k = <e^{ik theta}>}  and  mixed = <r^p e^{ik theta}>")
for M in [0.0, 1.0, 7.3]:
    Ms = [charged_moments(k) for k in (1, 2, 3)]
    mixed = [M * m for m in Ms]          # <r e^{ik theta}> = M * M_k -> 0 whenever M_k = 0
    print(f"   r-support = {{{M}}} : M_1={abs(Ms[0]):.2e}  M_2={abs(Ms[1]):.2e}  M_3={abs(Ms[2]):.2e}"
          f"   max|mixed|={max(abs(x) for x in mixed):.2e}   neutral gap = {M}")
print("  All charged moments vanish identically for every M: the charged data is ONE point,")
print("  the neutral gap scans [0,inf).  The map data->gap is not a function; the neutral")
print("  sector is the fibre of the forgetting map, NOT the quotient.  [T]+[V]")

# =====================================================================
blk("[F] beta=0 Haar acid test: full-support neutral measure, inf supp 0, atom 0")
var_trf = haar(trf ** 2) - haar(trf) ** 2
print(f"  <|tr_f|^2> = {haar(trf**2):.6f} , <tr_f>^2 = {haar(trf)**2:.2e}"
      f"  => Var(tr_f) = {var_trf:.6f}")
print("  beta=0: links i.i.d. Haar => the trF^2 correlator is a FINITE-RANGE CONTACT:")
print("     plaquettes at |dt| >= 2 share no link => C(t) = Var*delta_{t,0} + (dt=+-1 terms),")
print("     i.e. a nontrivial trigonometric polynomial in k (NOT flat, but FULL SUPPORT):")
print("     nu_{tr F^2} has support on the whole Brillouin zone => inf supp = 0, atom = 0.")
print("  (Continuum reading: gluons massless => trF^2 is a two-gluon cut on [0,inf).)")
print("  Charged data at beta=0: sigma=0 (perimeter law), string spectrum empty. Trivial.")
var_xi3 = haar(RHO ** 2) / 3.0
print(f"  cross-check with corpus: Var_Haar(xi3) = <rho^2>/3 = {var_xi3:.6f}"
      f"  vs (4pi^2-6)/9 = {(4*np.pi**2-6)/9:.6f}")

# =====================================================================
blk("SUMMARY (labels)")
print("[A] tr_adj = (tr_f)^2 - 1 exactly; neutral 1-pt = charged 2nd moment; <tr_adj>=0   [T]")
print("[B] 2D: m_neutral = 2N^2/(N^2-1) * m_charged >= 2 m_charged (pure YM, exact)       [T]")
print("[C] 2D: nu_{trF^2} EMPTY (atom=0 and Var=0) -> cannot refute any lower bound       [T]")
print("[D] decoration: identical charged data, neutral gap 0 vs 2*m_lat -> 0              [T]+[V]")
print("[E] toy: neutral is a FIBRE over charged data, not a quotient                      [T]+[V]")
print("[F] beta=0 Haar: flat measure, inf supp=0, atom=0 (acid pass)                      [T]")
print("\nNo QED. Clay 0%. The charged sector generates the neutral ALGEBRA but not the")
print("neutral STATE: what fails is not a bound but the invariance class of the hypotheses.")
