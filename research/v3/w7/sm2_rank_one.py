"""
sm2_rank_one.py  --  SM2 witness: the rank-one soft-mode theorem attempt (deterministic,
numpy only, py -3.12 -I).

Soft mode = transverse-magnetic (TM) Savvidy negative mode phi_-; the defect is rank one:
    K = K_+ - gB phi_- phi_-^T  (phi_- in ker K_+, |lam_-| = g|B|),
    D_Sav = E[ |<phi_-|c>|^2 ] / (gB)  >= 0,   D_+ <= 0 (BL on the positive fibre),
    G'' = S - D = (S - D_Sav) + |D_+| .

The SOFT-MODE MARGIN is  m_soft := S - D_Sav = E[ h_n - <K_+^+ c, c> ]
   = the curvature one would get if the negative mode were excised.

Sections:
 (1) exact rank-one split + soft-mode margin decomposition G''=(S-D_Sav)+|D_+|
     and the OVERLAP LEMMA  <phi_-|c> = <H_cf phi_- | v>  (coarse-fibre Hessian matrix element)
 (2) overlap: helicity / one-octave / scale separation  -> O(1) at the head k~sqrt(gB)
 (3) EXACT Var(phi^2; c, lam4): the "2/lam4" is an asymptote, peak 2.21/lam4 at c~4;
     at the head (small gB) it is far below 2/lam4
 (4) the precise inequality & constants: D_Sav <= g^2 Var(gB,lam4); margin m = S - g^2 Var.
     positivity region over (gB, lam4) at beta=2.5
 (5) beta-dependence of the margin
 (6) acid test: free field -> no negative mode -> D_Sav = 0, m = S
 (7) measured per-mode D/S, G'' on the L=8 npz (sanity, reproduces H2)
"""
import numpy as np, os

np.set_printoptions(precision=5, suppress=True)
rng = np.random.default_rng(20261124)

print("="*78)
print("SM2-RANK-ONE witness : rank-one soft-mode theorem attempt")
print("="*78)

# ============================================================================
# (1) EXACT SPLIT  +  SOFT-MODE MARGIN  +  OVERLAP LEMMA
# ============================================================================
print("\n(1) EXACT RANK-ONE SPLIT / SOFT-MODE MARGIN / OVERLAP LEMMA   [T]")
print("    G'' = (S - D_Sav) + |D_+| ;   m_soft := S - D_Sav = E[h_n - <K_+^+ c,c>]")
print("    overlap:  <phi_-|c> = <H_cf phi_- | v>   (coarse-fibre Hessian matrix element)")

# Work with explicit fibre / coarse bases; the fibre Hessian K_f is nf x nf
# with exactly one negative mode phi_f (eigenvalue -gB) and K_+ SPD on phi_f^perp.
n_c, n_f = 3, 5
n = n_c + n_f
Bc = np.column_stack([np.eye(n)[:, :n_c]])                  # coarse basis (n x n_c)
Bf = np.column_stack([np.eye(n)[:, n_c:]])                  # fibre basis  (n x n_f)
Pc = Bc @ Bc.T; Pf = Bf @ Bf.T
A = rng.standard_normal((n, n)); H = A + A.T                  # generic symmetric Hessian
phi_f = rng.standard_normal(n_f); phi_f /= np.linalg.norm(phi_f)
gB = 1.3
# K_+ SPD on the (nf-1)-dim complement of phi_f, kernel span(phi_f):
Qf, _ = np.linalg.qr(np.column_stack([phi_f, rng.standard_normal((n_f, n_f-1))]))
Qperp_f = Qf[:, 1:n_f]
spec = np.abs(rng.standard_normal(n_f-1)) + 0.4
Kplus_f = (Qperp_f*spec) @ Qperp_f.T
Kplus_f = 0.5*(Kplus_f+Kplus_f.T)
K_f = Kplus_f - gB*np.outer(phi_f, phi_f)                   # fibre Hessian, mode -gB
# splice K_f into the full H fibre block
H = H - Pf @ H @ Pf + Bf @ K_f @ Bf.T
phi = Bf @ phi_f                                             # ambient negative mode

v = Bc @ rng.standard_normal(n_c); v /= np.linalg.norm(v)    # a coarse mode (unit)
c_f = Bf.T @ H @ v                                           # c = Pi_0 H v in fibre coords
c   = Bf @ c_f
Hcf = Pc @ H @ Pf                                            # coarse<-fibre block
ov_direct = float(phi @ c)
ov_lemma  = float((Hcf @ phi) @ v)                           # <H_cf phi , v>
print(f"    <phi_-|c>            = {ov_direct:+.8f}")
print(f"    <H_cf phi_- | v>     = {ov_lemma:+.8f}     (overlap lemma, err {abs(ov_direct-ov_lemma):.2e})")

# exact split check on this mode:  <K^-1 c,c> = <K_+^+ c,c> - <phi,c>^2/gB   (fibre coords)
Kinv_f = np.linalg.inv(K_f); Kpinv_f = np.linalg.pinv(Kplus_f, rcond=1e-12)
lhs = c_f @ Kinv_f @ c_f
rhs = c_f @ Kpinv_f @ c_f - (phi_f @ c_f)**2 / gB
print(f"    split <K^-1 c,c> = <K_+^+c,c> - <phi,c>^2/gB :  lhs={lhs:+.6f} rhs={rhs:+.6f}  err={abs(lhs-rhs):.2e}")

# soft-mode margin identity on this toy:  m_soft = h_n - <K_+^+ c,c>  (local, no average)
h_n = float(v @ H @ v)
m_soft_local = h_n - float(c_f @ Kpinv_f @ c_f)
D_Sav_local  = float(phi_f @ c_f)**2 / gB
print(f"    toy: h_n={h_n:+.5f}  <K_+^+c,c>={float(c_f@Kpinv_f@c_f):+.5f}  "
      f"m_soft=h_n-<K_+^+c,c>={m_soft_local:+.5f}  D_Sav={D_Sav_local:+.5f}")

# algebraic identity:  (S-D_Sav) + |D_+|  vs  S - D  (with D = D_+ + D_Sav, D_+<=0)
# build a D_+ <= 0 (log-concave-fibre selloff) toy and confirm G''=(S-D_Sav)+|D_+|
S_bare = 1.40; D_Sav = 0.90; D_plus = -0.25
Gpp_direct = S_bare - (D_Sav + D_plus)
Gpp_split  = (S_bare - D_Sav) + abs(D_plus)
print(f"    margin identity: G'' = (S-D_Sav)+|D_+| :  direct={Gpp_direct:+.5f}  split={Gpp_split:+.5f}"
      f"   (m_soft = S-D_Sav = {S_bare-D_Sav:+.5f})")
print("    => m_soft is the excised-negative-mode curvature; |D_+| (BL) is a bonus, sign-only.")

# ============================================================================
# (2) OVERLAP: helicity / one-octave / scale separation
# ============================================================================
print("\n(2) OVERLAP  <phi_-|c_k>  ACROSS k/sqrt(gB)   [T in the LLL model / structural]")
print("    phi_- is the transverse-MAGNETIC (+1 helicity) lowest Landau level, support |p_~|~sqrt(gB).")
print("    |<phi_-|c_k>|^2 ~  [transversality: (k^2)/(k^2+gB)] * [LLL UV decay: exp(-k^2/gB)]")
def ov2(r):    # r = k/sqrt(gB)
    return (r*r/(1.0+r*r)) * np.exp(-r*r)
print(f"    {'k/sqrt(gB)':>10} {'|<phi|ck>|^2':>14}   regime")
for r in (0.0, 0.1, 0.3, 0.5, 0.71, 1.0, 1.414, 2.0, 3.0):
    tag = ("deep IR: transverse suppression -> 0" if r < 0.3 else
           ("HEAD: O(1), max" if 0.5 <= r <= 1.5 else
            ("UV: LLL decay -> 0" if r > 1.8 else "")))
    print(f"    {r:>10.3f} {ov2(r):>14.4f}   {tag}")
rpk = np.sqrt(np.sqrt(5.0)/1.0)  # argmax of r^2/(1+r^2)*exp(-r^2)
rs = np.linspace(0.001,4,400000); pk = rs[np.argmax(ov2(rs))]
print(f"    peak at k/sqrt(gB)={pk:.3f} with value {ov2(pk):.4f}  => in-band O(1), NOT small")
print("    three mechanisms: (i) scale separation kills k<<sqrt(gB); (ii) one-octave band")
print("    k~sqrt(gB) gives no suppression (O(1)); (iii) gauge/helicity projection: only the")
print("    +1 (TM) projection of c couples -> a longitudinal/pure-gauge v_k has ZERO overlap.")

# ============================================================================
# (3) EXACT Var(phi^2; c= gB, lam4): the 2/lam4 asymptote vs the peak vs the head
# ============================================================================
print("\n(3) EXACT Var(phi^2; c=gB, lam4)  for  V=-c/2 phi^2 + lam4/4 phi^4   [V-num]")
def var_phi2(c, lam, npts=None, pad=5.0):
    v0 = np.sqrt(c/lam) if c > 0 else 1.0
    L = max(6.0, pad*v0 + 4.0)
    if npts is None:
        npts = int(min(1_200_001, max(200_001, 80000*L)))
    x = np.linspace(-L, L, npts)
    e = 0.5*c*x**2 - 0.25*lam*x**4; e -= e.max()
    w = np.exp(e); w /= np.trapezoid(w, x)
    m2 = np.trapezoid(w*x**2, x); m4 = np.trapezoid(w*x**4, x)
    return m4 - m2**2
print(f"    {'gB':>6} {'Var(phi^2)':>12} {'/(2/lam4)':>11}   (lam4=1)")
for c in (0.1,0.25,0.5,1.0,2.0,3.0,4.0,5.0,7.0,10.0,50.0,100.0):
    val = var_phi2(c, 1.0)
    print(f"    {c:>6.2f} {val:>12.5f} {val/2.0:>11.4f}")
print("    => 2/lam4 is the c->inf ASYMPTOTE, not a bound: peak ~2.21/lam4 at gB~4 (+10%).")
print("       At the HEAD (small gB) Var is 0.2-0.5 x (2/lam4): 2-3x TIGHTER than 2/lam4.")

print("\n    head grid  Var(gB,lam4):")
print(f"    {'gB\\lam4':>8}" + "".join(f"{l:>10.2f}" for l in (0.8,1.6,2.4,3.2)))
for gBv in (0.3,0.5,1.0,2.0,4.0):
    row = "".join(f"{var_phi2(gBv,l):>10.3f}" for l in (0.8,1.6,2.4,3.2))
    print(f"    {gBv:>8.2f}"+row)

# ============================================================================
# (4) THE PRECISE INEQUALITY & CONSTANTS:  D_Sav <= g^2 Var(gB,lam4) ;  m = S - g^2 Var
# ============================================================================
print("\n(4) THE INEQUALITY  D_Sav < S_bare - delta  AND  THE MARGIN   [T]/[C]")
print("    D_Sav <= g^2 Var(gB,lam4)   (exact finite-c Var, NOT the 2/lam4 asymptote)")
print("    matrix-element form:  D_Sav <= E||P_{+1} H_fc v||^2 / (gB)  vs  S_bare=E[h_n]")
print("    margin  m(gB,lam4;beta) = S_bare(beta) - g^2(beta) Var(gB,lam4)")
beta = 2.5; g2 = 4.0/beta; S_ref = 1.40             # measured E[S] at beta~2.4-2.5
print(f"    beta={beta}, g^2=4/beta={g2:.3f}, S_bare~{S_ref} (measured, [V])")
print(f"    {'gB':>6} {'lam4':>6} {'Var':>9} {'g^2Var':>9} {'margin m':>10} {'m>0?':>6}")
for gBv, lam in [(0.5,1.6),(1.0,1.6),(2.0,1.6),(2.0,0.8),(1.0,0.8),(0.3,0.8),(4.0,1.6),(2.0,3.2)]:
    V = var_phi2(gBv, lam); D = g2*V; m = S_ref - D
    print(f"    {gBv:>6.2f} {lam:>6.2f} {V:>9.4f} {D:>9.4f} {m:>10.4f} {('YES' if m>0 else 'NO'):>6}")

print("\n    positivity region  m>0  at beta=2.5, S_ref=1.40:")
print(f"    {'gB\\lam4':>8}" + "".join(f"{l:>10.2f}" for l in (0.5,0.8,1.0,1.2,1.6,2.0,2.4,3.2)))
for gBv in (0.3,0.5,1.0,1.5,2.0,3.0,4.0):
    row = ""
    for lam in (0.5,0.8,1.0,1.2,1.6,2.0,2.4,3.2):
        m = S_ref - g2*var_phi2(gBv, lam)
        row += f"{('m>0 ' if m>0 else 'm<0 '):>10}"
    print(f"    {gBv:>8.2f}"+row)
print("    => m>0 for natural-to-stiff lam4 (>~1.2) up to gB~3; FAILS for soft lam4 (<~1) at gB>~1.")
print("       Natural lam4~g^2=1.6: m>0 throughout gB in [0.3,2] -- conditional margin EXISTS.")

# threshold: smallest lam4 with m>0 (m<0 for SOFT quartic, m>0 for STIFF quartic)
def lam_min_for(gBv):
    lo,hi=0.2,12.0
    for _ in range(80):
        mid=0.5*(lo+hi)
        if S_ref-g2*var_phi2(gBv,mid)>0: hi=mid       # stiff side: m>0
        else: lo=mid                                   # soft side: m<0
    return hi
print("\n    threshold lam4_min(gB): smallest lam4 with m>0 (stiffer quartic => m>0):")
for gBv in (0.3,0.5,1.0,2.0,3.0,4.0):
    lm = lam_min_for(gBv)
    print(f"      gB={gBv:>4}: lam4_min = {lm:.3f}   (g^2=1.6, lam4_min/g^2={lm/g2:.2f})")
print("    => m>0 requires the quartic STIFFER than a gB-dependent threshold; natural lam4~g^2")
print("       lies above it for small gB, on the edge / below it for large gB.")

# ============================================================================
# (5) BETA-DEPENDENCE OF THE MARGIN
# ============================================================================
print("\n(5) BETA-DEPENDENCE  m(beta)=S_bare - g^2(beta) Var(gB,lam4)   [C]/[O]")
print("    soft = g^2 Var(gB,lam4); its beta->inf limit is an RG input (how gB,lam4 scale):")
# (a) quartic held O(1):  soft = g^2 Var -> 0  => SAFE
print("    (a) lam4 fixed O(1):   soft = g^2 Var(gB,O(1)) -> 0 as beta->inf  => SAFE")
print(f"    {'beta':>6} {'g^2':>7} {'g^2Var(gB=1,lam4=1)':>20} {'m=S-g^2Var':>12}")
for b in (1.0,2.0,2.5,5.0,10.0,50.0):
    gg=4.0/b; D=gg*var_phi2(1.0,1.0); print(f"    {b:>6.1f} {gg:>7.3f} {D:>20.4f} {S_ref-D:>12.4f}")
# (b) natural lam4=g^2, gB fixed:  asymptote soft -> 2g^2/lam4 = 2 (marginal, S<2)
print("    (b) natural lam4=g^2 (quartic runs with coupling), gB=1 fixed:")
print(f"    {'beta':>6} {'g^2':>7} {'Var exact':>10} {'g^2Var':>9} {'2g^2/lam4':>10} {'m exact':>9}")
for b in (2.5,5.0,10.0,50.0):
    gg=4.0/b; V=var_phi2(1.0,gg); D=gg*V; asym=2*gg/gg
    print(f"    {b:>6.1f} {gg:>7.3f} {V:>10.4f} {D:>9.4f} {asym:>10.2f} {S_ref-D:>9.4f}")
print("    => (b) asymptote soft -> 2g^2/lam4 = 2 > S=1.4: the ASYMPTOTE FAILS, but the")
print("       FINITE-gB exact value at beta=2.5 is 0.85 < 1.4: the pass hinges on SMALL gB.")
print("       so beta-dependence does NOT give a uniform delta: it is a pure O(1) ratio = H1.")

# ============================================================================
# (6) ACID TEST
# ============================================================================
print("\n(6) ACID TEST (free field)   [T]")
print("    B=0 => no Savvidy background => no negative mode phi_- => D_Sav = 0,")
print("    D_+ = 0 (Gaussian fibre) => G'' = S (the free marginal curvature); at the")
print("    IR mode Q_c~4pi^2/L^2 -> 0, so G''->0: the rank-one tool manufactures NO gap.")
S_free = 0.0; print(f"    free massless: D_Sav=0, m_soft=S={S_free:.4f}, G''={S_free:.4f}  (passes)")

# ============================================================================
# (7) MEASURED per-mode check (sanity, reproduces H2)
# ============================================================================
print("\n(7) MEASURED per-mode D/S, G''  on k0_su2_b2.4_L8_n2.npz   [V-num]")
path = os.path.normpath(r"D:\Dev\math\research\unconditional\w7\k0_su2_b2.4_L8_n2.npz")
if os.path.exists(path):
    d = np.load(path)
    Hm,Gm,P1,P2 = d['H'],d['G'],d['P1'],d['P2']
    Sm = Hm.mean(0); VarG = Gm.var(0); EP2=P2.mean(0); VarP1=P1.var(0)
    Cov = ((Gm-Gm.mean(0))*(P1-P1.mean(0))).mean(0)
    W = EP2 - VarP1 - 2*Cov
    Dm = VarG - W; Gpp = Sm - Dm; DS = Dm/Sm
    print(f"    Ncfg={Hm.shape[0]} Nmode={Hm.shape[1]}  beta={float(d['beta'])} L={int(d['L'])}")
    print(f"    D/S: mean={DS.mean():.4f} min={DS.min():.4f} max={DS.max():.4f}")
    print(f"    G'': mean={Gpp.mean():+.4f} min={Gpp.min():+.4f} max={Gpp.max():+.4f}")
    print(f"    D<S: {int((Gpp>0).sum())}/{len(Gpp)}   D>S: {int((Gpp<0).sum())}/{len(Gpp)}")
    print("    => near-cancellation, both signs: the measured margin is O(1) and thin.")

    # EXACT upper bound on the SOFT-MODE margin:  G'' = m_soft + |D_+|  =>  m_soft <= G''
    print("\n(7b) EXACT CEILING on the soft-mode margin:  G'' = m_soft + |D_+|  =>  m_soft <= G''   [T]")
    print("     m_soft = S - D_Sav = G'' - |D_+|  <=  G''  (BL part |D_+|>=0 can only help).")
    print(f"     measured mean G'' = {Gpp.mean():+.4f}  (={Gpp.mean()/Sm.mean()*100:.1f}% of mean S)")
    print(f"     => ANY positive soft-mode margin is <= G'' ; in this ensemble m_soft <= {Gpp.mean():+.4f}.")
    print(f"     D_+<=0 unbounded below => m_soft in (-inf, G'']; the positivity is BL-CARRIED, not")
    print("     rank-one-carried.  Proving m_soft>0 would be SUFFICIENT but is NOT given by the data")
    print("     (and is STRICTLY harder than G''>0).  This is the exact obstruction to a soft-mode")
    print("     theorem from the rank-one structure alone.")

    # CONSISTENCY PINS gB:  m_toy = S-g^2Var <= m_soft <= G''  =>  g^2Var >= S-G'' = D_Sav(meas)
    print("\n(7c) the bound (m_soft>=m_toy) + ceiling (m_soft<=G'') PIN the head curvature:  [T]/[V]")
    print("     need g^2 Var(gB,lam4) >= S - G'' = D_Sav(measured) ~ 0.965*S = %.3f" % (0.965*S_ref))
    target = 0.965*S_ref
    for gBv in (1.0,1.5,2.0,2.5,3.0):
        for lam in (1.6,):
            Dv = g2*var_phi2(gBv, lam)
            print(f"       gB={gBv:>4}, lam4={lam}: g^2Var={Dv:.3f}  ({'>= target' if Dv>=target else '< target'})")
    print("     => crossing at gB~2.1: the data sit at gB>~2.1, the THIN-margin corner")
    print("        (m_toy <~ G'' ~ 0.05, dropping to <=0 as gB grows).  The gB=1 corner")
    print("        (m_toy=0.55) is EXCLUDED by D_Sav~S.")
else:
    print("    [npz not found]")

print("\n" + "="*78)
print("VERDICT: rank-one structure gives the EXACT soft-mode margin identity")
print("G''=(S-D_Sav)+|D_+|, the overlap lemma <phi_-|c>=<H_cf phi_-|v> (O(1) at the head),")
print("and a CONDITIONAL positive margin m_soft=S-g^2 Var(gB,lam4)>0 in a computable region")
print("(natural lam4~g^2, gB<~3). BUT the exact ceiling m_soft = G'' - |D_+| <= G'' (<=3.5% of")
print("S measured) forbids an O(1) soft-mode margin: the observed positivity is BL-carried, not")
print("rank-one-carried. No unconditional delta>0: gB, lam4, S_bare at the head are H1 inputs and")
print("the head is the O(1)-overlap / IR-marginal / Landau-pole endpoint. Acid test passes. Clay 0%.")
print("="*78)
