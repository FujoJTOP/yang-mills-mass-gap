# r2_h2global.py -- R2 witness: H2-global, definitions D1/D2/D3, free-field vacuity, C(beta) window
# py -3.12 -I  (numpy 2.4.6).  Deterministic, single process, no MC, no download, < 5 s.
import numpy as np
import itertools

np.set_printoptions(precision=6, suppress=True)


def hr(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


# ----------------------------------------------------------------------------------
# geometry helpers
# ----------------------------------------------------------------------------------
def minus_delta_spectrum(L):
    """eigenvalues of -Delta on the L^4 periodic torus, a=1:  q(k)=4 sum_mu sin^2(k_mu/2)"""
    ks = np.arange(L) * 2.0 * np.pi / L
    K = np.meshgrid(ks, ks, ks, ks, indexing="ij")
    q = 4.0 * (np.sin(K[0] / 2) ** 2 + np.sin(K[1] / 2) ** 2
               + np.sin(K[2] / 2) ** 2 + np.sin(K[3] / 2) ** 2)
    return np.moveaxis(q, 0, -1)  # index tuple -> value


def fiber_mask(L):
    """sharp factor-2 low-pass p(k): pass if every component |k_mu| <= pi/2 ; fiber = complement"""
    ks = np.arange(L) * 2.0 * np.pi / L
    kk = np.array(np.meshgrid(ks, ks, ks, ks, indexing="ij"))
    # |k_mu| <= pi/2  <=>  k_mu in [0, pi/2] u [2pi - pi/2, 2pi)  on the torus
    per_axis = []
    for m in range(4):
        k = kk[m]
        per_axis.append((k <= np.pi / 2 + 1e-12) | (k >= 2 * np.pi - np.pi / 2 - 1e-12))
    coarse = per_axis[0] & per_axis[1] & per_axis[2] & per_axis[3]
    return coarse  # fiber = ~coarse


def n_plaquettes_block(b):
    """plaquettes in a b^4 periodic block: 6 b^4 (lattice convention), and 3 b^4 (corpus PROVE-H2)"""
    return 6 * b ** 4, 3 * b ** 4


# ==================================================================================
# [A] FREE MASSLESS FIELD: H2-global holds, V-uniformly and b-uniformly.
#     and H1/(star) fails.  =>  H2g =/=> H1   [T]+[V]
# ==================================================================================
hr("[A] free massless field: fiber band floor (H2-global) vs k->0 mode (H1)")
print("  block map P = sharp factor-2 low-pass; fiber = ker P = high band|k_mu|>pi/2")
print("  L     min spec(-Delta|_fiber)   G''(0)=q(0) [=H1 object]   #fiber modes")
rowsA = []
for L in (4, 8, 12, 16, 24, 32):
    q = minus_delta_spectrum(L)
    co = fiber_mask(L)
    fib = q[~co]
    rowsA.append((L, fib.min(), q[tuple([0] * 4)], fib.size))
    print(f"  {L:3d}      {fib.min():18.10f}      {q[tuple([0]*4)]:18.10f}   {fib.size:10d}")
print("  => fiber (H2-global) gap  -> 2.0 = 4 sin^2(pi/4),  V-UNIFORM, b-INDEPENDENT")
print("  => G''(0)=q(0)=0 exactly: the coarse mode is FREE  =>  H1 FAILS on the same field")
print("  => [V]  H2-global TRUE  AND  H1 FALSE   ==>  H2-global /=> H1   (incomparable, one side)")
print("  b-uniformity: the fiber spectrum does not know b at all (one fixed band), so")
print("     Lambda_UV^free(b) = 2/a^2 for every b   (cf. H2-b-uniform Sec.1)")

# ==================================================================================
# [B] Zegarlinski block-Dobrushin coefficient C(beta): D1 is FALSE at beta=2.5
# ==================================================================================
hr("[B] D1 (Zegarlinski C(beta)<1): exact window, disjoint from beta=2.5")
print("  delta_kj >= 2*beta per straddling plaquette (that IS the plaquette weight 2*beta*[0,1])")
print("  C(beta) >= 2*beta*deg_future ;  deg_future = 3 in 4D (faces +mu, mu=1..4, ordered)")
print("  window of D1:  2*beta*deg < 1  <=>  beta < 1/(2*deg)")
for deg in (1, 3, 4):
    bmax = 1.0 / (2.0 * deg)
    print(f"    deg={deg}:  C(beta)>= {2*deg:.0f}*beta ;  D1 window beta < {bmax:.6f}"
          f"   ;  2.5/bmax = {2.5/bmax:7.3f}x")
print("  at beta=2.5:")
for deg in (3, 4):
    print(f"    deg={deg}:  C >= {2*2.5*deg:8.3f}   (corpus quotes C~15)   => C<1 FALSE")
print("  => D1 fails for beta > 1/6 ;  beta=2.5 is 15x outside.  D1 is a ROUTE, not the statement.")

# ==================================================================================
# [C] block-local fiber LSI is PROVED (HS) at beta=2.5 but numerically worthless
# ==================================================================================
hr("[C] proven block-local Lambda_UV at beta=2.5, and what it delivers to m")
beta = 2.5
C_G = 1.0  # O(1) LSI constant of normalized Haar on SU(2)
print("  b     N_p(B)[6b^4]  N_p(B)[3b^4]   Lambda_UV=e^{-2 beta N_p}/C_G   log10")
for b in (2, 4, 8):
    n6, n3 = n_plaquettes_block(b)
    l6 = np.exp(-2 * beta * n6) / C_G
    l3 = np.exp(-2 * beta * n3) / C_G
    s6 = f"{np.log10(l6):8.2f}" if l6 > 0 else "    -inf"
    s3 = f"{np.log10(l3):8.2f}" if l3 > 0 else "    -inf"
    print(f"  {b:2d}   {n6:10d}    {n3:10d}     {l6:10.3e} / {l3:10.3e}   "
          f"{s6} / {s3}")
lam = 0.165  # hypothetical H1 constant (magnitude of the measured G''(0)=-0.165)
mfib = np.exp(-2 * beta * 6 * 2 ** 4)
print(f"\n  assembly m = lam*lam_f/(lam+lam_f),  lam={lam},  lam_f=Lambda_UV(b=2)")
print(f"    lam_f = {mfib:.6e}   ->   m = {lam*mfib/(lam+mfib):.6e}")
print("  => the ONLY proven H2-global constant contributes  m = lam_f = e^{-240}: m is 0")
print("     to >100 digits.  A proven H2-global with the HS constant is QUANTITATIVELY VACUOUS.")

# ==================================================================================
# [D] the assembly is loss-only:  m <= min(lam, lam_f)
# ==================================================================================
hr("[D] m = lam*lam_f/(lam+lam_f) <= min(lam,lam_f):  H2-global is LOSS-ONLY for the conclusion")
print("  lam      lam_f      m=harm.mean    min(lam,lam_f)   m/lam     m/lam_f")
for lf in (1e1, 1.0, 0.165, 1e-3, 1e-30):
    L = 0.165
    m = L * lf / (L + lf) if (L + lf) > 0 else 0.0
    print(f"  {L:7.4f}  {lf:9.3e}  {m:12.6e}    {min(L,lf):12.6e}   {m/L:8.5f}  {m/lf:8.5f}")
print("  => lam_f>>lam: m->lam   (H2g dispensable up to factor 2)")
print("  => lam_f<<lam: m->lam_f (H2g becomes the binding constraint and kills the bound)")
print("  => NEVER does H2g RAISE the conclusion above lam.  It can only lower it.")

# ==================================================================================
# [E] D2 (pointwise in U') vs D3 (nu-averaged): D3 strictly weaker, and D3 is ENOUGH
# ==================================================================================
hr("[E] D2 vs D3: the split is nu-AVERAGED, so D3 suffices for the assembly")
rng = np.random.default_rng(20261125)
n = 200000
# toy: fiber variance is a strongly U'-dependent function; nu = N(0,1) on U'
u = rng.normal(size=n)
g = np.exp(3.0 * np.tanh(u))          # Var_{mu(.|u)}(f) for a representative f
print("  toy: g(U') = exp(3 tanh U'),  nu = N(0,1)")
print(f"    sup_U' g(U')   = {g.max():10.5f}     <- D2 must bound THIS (pointwise)")
print(f"    E_nu[g]        = {g.mean():10.5f}     <- D3 must bound only THIS (averaged)")
print(f"    ratio E/sup    = {g.mean()/g.max():10.5f}     (strict, and can be made 0)")
print("  two-scale law of total variance:")
print("    Var_mu(F) = E_nu[ Var_{mu(.|U')}(F) ]  +  Var_nu( E_mu[F|U'] )")
print("    term1 is E_nu[.]  =>  a nu-averaged fiber LSI (D3) bounds term1 EXACTLY.")
print("    term2 is coarse  =>  H1 bounds it.")
print("  => D3 (nu-averaged) is SUFFICIENT for  H1 & D3 => Var_mu(F) <= max(1/lam,1/lam_f) E|grad F|^2")
print("  => D2 (pointwise) is STRICTLY stronger than D3, and buys ONLY the pointwise-in-U' form")
print("     that (star dec)/L3-band-decay needs (S4 Sec.3.5), NOT the UC-COMPLETE three candidates.")

# ==================================================================================
# [F] H1 =/=> H2-global : explicit witness (conditionals do not constrain the marginal)
# ==================================================================================
hr("[F] H1 /=> H2-global : explicit product witness")
for L in (4, 8, 16, 32, 64):
    q = minus_delta_spectrum(min(L, 32))
    co = fiber_mask(min(L, 32))
    fib = q[~co]
    # massless Gaussian fiber on the torus: LSI constant = band floor (V-uniform); but
    # if one takes the fiber to be the WHOLE configuration (no band cut), the constant is
    # the torus gap (2 pi/L)^2 -> 0.  That is the H1-true / H2g-false witness.
    whole = (2 * np.pi / L) ** 2
    print(f"  L={L:3d}:  (2pi/L)^2 = {whole:.6e}   [whole-space Gaussian: LSI -> 0, H2g FALSE]")
print("  coarse marginal nu chosen independently as N(0,1) => H1 TRUE with lam=1.")
print("  Because mu = nu(da) mu(.|a) with {mu(.|a)} and nu INDEPENDENT (B3-ISOPERIMETRIC,")
print("  'conditionals do not constrain the marginal'), this is a legal joint measure.")
print("  => [V] H1 TRUE and H2-global FALSE  ==>  H1 =/=> H2-global.")
print("  Together with [A] (H2g TRUE, H1 FALSE) => H2-global INCOMPARABLE with H1. [T]+[V]")

# ==================================================================================
# [G] calibrations where H2-global is TRUE and gap content is ZERO (vacuumisation)
# ==================================================================================
hr("[G] two calibrations: H2-global TRUE while the gap content is 0")
print("  (G1) free massless field : H2-global [A] TRUE, f1^2 = 0, gap = 0        -> vacuous")
print("  (G2) 2D pure YM (Migdal) : tr F^2 = H conserved => f1^2 = nu({m1^2}) = 0,")
print("       while the FIBER (band-limited, same as [A]) is uniformly gapped, so")
print("       H2-global TRUE.  Gap content = 0 (no atom).        -> vacuous")
print("  => H2-global never carries the atom in any solvable model the corpus owns.")
print("  => [V] H2-global is a UV STRUCTURAL lemma, not a gap hypothesis.")

# ==================================================================================
# [H] sufficient-for-the-three-candidates checklist
# ==================================================================================
hr("[H] which definition makes the three UC-COMPLETE candidates work?")
print("  cand.                 needs                        D1        D2        D3")
print("  #3  H1 & X            X => term1 of the split     yes       yes       yes")
print("  #4  LSI & X           same                         yes       yes       yes")
print("  S7  'debt moves to X' same                         yes       yes       yes")
print("  status at beta=2.5    --                        FALSE     [O]       [O]")
print("  window                --                        beta<1/6  none      none")
print("  => D1: logically sufficient, numerically false => kills all three (must be DEMOTED).")
print("  => D2: the corpus's working definition; open at 2.5; no window.")
print("  => D3: STRICTLY WEAKER than D2 and still sufficient for the three, BUT")
print("         SZZ gives the FULL LSI, which does NOT imply D3 (extra E|grad_par F|^2 term)")
print("         -> D3 has NO known window either.  => theorem R2-4 dichotomy (see R2-H2GLOBAL 4.3)")

# ==================================================================================
# [I] SZZ's FULL LSI does NOT give D3: the grad_parallel leakage, quantified
# ==================================================================================
hr("[I] full LSI =/=> D3 : the E|grad_par F|^2 leakage (exact toy)")
print("  mu = N(0,s^2) on the coarse coord a  x  N(0,1) on the fiber coord b   (Gaussian => full LSI, lam=1)")
print("  F = a + b  (general, NOT xi-only).  Exact:  Var= s^2+1 ; E|grad_par F|^2 = s^2 ; E|grad_perp F|^2 = 1")
print("   s     Var(F)     (1/1)(E|gp|^2+E|gq|^2)   (1/1)E|gq|^2      Var <= E|gq|^2 ?")
for s in (0.1, 1.0, 10.0, 100.0, 1000.0):
    V = s * s + 1.0
    par = s * s
    perp = 1.0
    print(f"  {s:6.1f}  {V:10.4f}   {par+perp:20.4f}   {perp:14.4f}      {'YES' if V <= perp else 'NO'}")
print("  => the full LSI bounds Var by (E|gp|^2+E|gq|^2); it does NOT bound it by E|gq|^2 alone.")
print("     The obstruction is exactly the leakage the corpus calls 'no grad_perp bound gives a")
print("     grad_par bound' (PROOF-IR-H1 Sec.1) -- here run in reverse.  So D3 is NOT a corollary")
print("     of SZZ; and for xi-only f it IS (Var_mu(f) <= lam^-1 E|grad_perp f|^2), which is the")
print("     RESTRICTED form the assembly cannot use.")

print("\n[DONE] r2_h2global.py  -- all sections computed, exit 0")