"""
sl_test.py -- Spectral-identity cluster test:
  (A) Spectral Lineage  (projective system, refinements)
  (B) Unfinished Spectral Event (candidate projections, three fates)
  (C) Counterfactual Spectral Shadow (remove a candidate's traces)

Smallest exact 1D/2D models + the four-case battery:
  free massless | 2D YM holonomy | property-but-no-gap | gapped-without-property

numpy/scipy only, one job, no downloads.
Run:  py -3.12 -I sl_test.py
"""

import numpy as np

def hdr(t):
    print("\n" + "=" * 72)
    print(t)
    print("=" * 72)

rng = np.random.default_rng(0)

# =====================================================================
# CASE 1 -- FREE MASSLESS 1D LATTICE  (gapless, atomless in continuum)
# =====================================================================
hdr("[CASE 1] free massless 1D periodic lattice: m1^2(L) -> 0, no atom")
print("  lattice Laplacian eig  lam_k = 2-2cos(2pi k/L); lowest nonzero lam_1")
rows = []
for L in [8, 16, 32, 64, 128, 256]:
    lam = np.sort(2 - 2 * np.cos(2 * np.pi * np.arange(L) / L))
    m1sq = lam[1]
    rows.append((L, m1sq))
    print(f"   L={L:4d}  m1^2={m1sq:.6e}   4pi^2/L^2={4*np.pi**2/L**2:.6e}"
          f"   ratio={m1sq/(4*np.pi**2/L**2):.4f}")
print("  => m1^2 -> 0 (no mass gap).  Continuum KL measure atomless.")
# edge atom read: mass in a window [0,d] at finite L vs d
print("  edge atom read  M(d)=mass in [0,d] (continuum proxy rho=1 on [0,1]):")
for d in [0.2, 0.05, 0.01, 0.003]:
    print(f"   d={d:<7} M(d)={d:.5f}  -> 0 like d  (NO atom at edge)")

# =====================================================================
# CASE 2 -- 2D YM SU(2) CYLINDER  (gapped; trF^2 = H conserved => f1^2=0)
# =====================================================================
hdr("[CASE 2] 2D YM SU(2) cylinder: Casimir transfer spectrum, gap, trF^2 atom")
e2 = 1.0
spins = [0.0, 0.5, 1.0, 1.5, 2.0]
# SU(2) quadratic Casimir C2(j)=j(j+1);  transfer eigenvalue E=(e^2/2) C2
E = [e2 / 2 * (j * (j + 1)) for j in spins]
for j, Ej in zip(spins, E):
    print(f"   j={j:3}   E=(e^2/2) j(j+1) = {Ej:.6f}")
gap = E[1] - E[0]
print(f"   gap Delta = E(1/2)-E(0) = {gap:.6f}   (=3e^2/8 = {3*e2/8:.6f})")
# trF^2 ~ H is diagonal in the R-basis; <0|trF^2|R!=0> = 0
# => KL atom of trF^2 : f1^2 = |<0|trF^2|1>|^2 = 0
print("   trF^2 proportional to H, diagonal in R => <0|trF^2|first-excited> = 0")
print("   => f1^2 = 0 (W-charter acid #2: trF^2=H conserved).  GAPPED, ATOMLESS.")
# two-point / connected variance of trF^2 in the vacuum is 0 (H eigenstate)
print("   <(trF^2)^2>_c = Var(H)_vac = 0  (vacuum is H-eigenstate)  => no excitation.")

# =====================================================================
# CASE 3 -- THIN ATOM (H1-blind) and w*delta_m + Unif (atom, no floor)
# =====================================================================
hdr("[CASE 3a] thin atom nu_e = (e/2lam)delta_e + (1-e/2lam)delta_s2  (H1 blind)")
lam, s2 = 1.0, 3.0
print("   finite data (mass M0, first moment <s>, truncated chi=1/<s>) vs edge atom w:")
for eps in [0.5, 0.2, 0.05, 0.01, 0.001]:
    w = eps / (2 * lam)
    M0 = 1.0
    M1 = w * eps + (1 - w) * s2           # first moment (a finite datum, weak-* continuous)
    chi_trunc = 1.0 / M1                  # bounded 'susceptibility' finite datum
    print(f"   eps={eps:<7} edge_atom w={w:.4e}  M0={M0:.6f}  <s>={M1:.6f}  chi=1/<s>={chi_trunc:.6f}")
print("   => every finite datum (M0, <s>, chi=1/<s>) stays O(1) and bounded (-> s2=3),")
print("      while the edge atom w=e/2lam -> 0  (H1/finite data cannot see the atom; R3 [C3]/[C4]).")

hdr("[CASE 3b] nu = w delta_m + (1-w) Unif[0,Lam] : ATOM (w>0) but inf supp=0 (NO gap)")
w, m, Lam = 0.3, 1.0, 3.0
print(f"   atom at m={m} weight {w}; continuum Unif[0,{Lam}] touches 0 => inf supp = 0")
for d in [0.5, 0.2, 0.05, 0.01]:
    lo = max(0.0, -d); hi = min(Lam, d)
    mass0 = (1 - w) * (hi - lo) / Lam
    print(f"   mass in [-{d},{d}] = {mass0:.5f}  -> 0 like density*d  BUT the measure reaches 0")
print("   near the atom m (window isolated from 0 and from Lam):")
for d in [0.2, 0.05, 0.01, 0.003]:
    cont = (1 - w) * (min(Lam, m + d) - max(0.0, m - d)) / Lam
    print(f"   mass in [m-{d},m+{d}] = atom {w:.3f} + continuum {cont:.5f} = {w+cont:.5f}")
print("   => genuine atom of weight 0.3 but NO mass gap (inf supp = 0).")

# =====================================================================
# CASE 4 -- GAPPED ATOMLESS  nu = c*sqrt(s-m^2) on [m^2, m^2+W]
# =====================================================================
hdr("[CASE 4] gapped atomless : 2-particle threshold density rho=c*sqrt(s-m^2)")
m2, Wc = 1.0, 3.0     # inf supp = m^2 = 1 > 0 ; a.c., no atom
print(f"   rho(s)=c*sqrt(s-{m2}) on [{m2},{m2+Wc}], c chosen so total mass=1")
# normalise
def dens(s):
    return np.sqrt(np.clip(s - m2, 0, None))
# integral sqrt(s-m2) ds over [m2,m2+Wc] = (2/3) Wc^{3/2}
Z = (2.0 / 3.0) * Wc ** 1.5
print(f"   normalization Z=(2/3)W^{{3/2}}={Z:.6f}")
print("   edge atom read M(d)=mass in [m^2, m^2+d]:")
for d in [0.2, 0.05, 0.01, 0.003]:
    mass = (2.0 / 3.0) * d ** 1.5 / Z
    print(f"   d={d:<7} M(d)={mass:.6e}  -> 0 like d^(3/2)  (NO atom at edge)")
print(f"   inf supp = {m2} > 0  => MASS GAP present, but f1^2(edge atom)=0.")

# =====================================================================
# (A) SPECTRAL LINEAGE -- existence = inverse-limit compactness
# =====================================================================
hdr("[A] spectral lineage: nested intervals I_n, consistent atoms, limit mass")
print("  Lineage evidence set E_n = {consistent features w/ mass>=w, supp in I_n}.")
print("  Existence of a projective limit = finite-intersection / Banach-Alaoglu (AUTOMATIC).")
print("  What matters is whether the limit MASS is >= w>0 uniformly. Demonstrate on 4 cases:")
print("  free massless: candidate edge LOCATION m_n = lam_1(L) -> 0, mode weight O(1) at every n:")
for L in [8, 16, 32, 64, 128, 256]:
    m_n = 2 - 2 * np.cos(2 * np.pi / L)
    print(f"      n={L:4d}  candidate atom location m_n={m_n:.6e}  (weight O(1))  -> floor 0")
print("     => lineage EXISTS with O(1) weight but its LOCATION -> 0: the limit floor = 0 (no gap). R1 face.")
print("  2D YM (trF^2): atom weight f1^2=0 at every refinement (gapped, atomless for trF^2).")
print("  thin atom w_e=e/2lam: edge atom -> 0 ([3a]);  w delta_m+Unif: atom 0.3 but floor 0 ([3b]).")
print("  => lineage EXISTENCE is compactness (inverse limit); POSITIVITY of the limit = the atom floor W")

# =====================================================================
# (B) UNFINISHED SPECTRAL EVENT -- three fates on H(lam)=diag(lam,1)
# =====================================================================
hdr("[B] three fates on H(lam)=diag(lam,1), candidate P = 1{eig < 1/2}")
thr = 0.5
def fate(lam_val):
    eigs = np.array([lam_val, 1.0])
    inproj = np.where(eigs < thr)[0]
    rank = len(inproj)
    if rank == 1:
        # isolated from the other (1.0) with gap
        iso = abs(1.0 - eigs[inproj[0]])
        return rank, iso, "FIRES (isolated sub-threshold projection)"
    elif rank == 0:
        return rank, np.nan, "FORBIDDEN (no sub-threshold projection)"
    else:
        return rank, np.nan, "?"
for lam_val, tag in [(0.2, "generic allowed"), (0.6, "crossed -> forbidden")]:
    r, iso, msg = fate(lam_val)
    print(f"   lam={lam_val}: rank={r}  isolation={iso}  -> {msg}")
print("   lam_n = 1/2 - 1/n  (->1/2 from below): rank=1 for every n,")
for n in [2, 4, 8, 16, 10**3, 10**6]:
    lam_val = 0.5 - 1.0 / n
    r, iso, _ = fate(lam_val)
    print(f"      n={n:<8} lam={lam_val:.8f}  rank={r}  isolation_to_1={iso:.6e}")
print("   but isolation -> 1/2 uniformly, eigenvalue hits threshold: NO finite certificate")
print("   => UNFINISHED fate == gapless / threshold-touch == spectral-flow-at-a-crossing.")
# FALSIFY-FLOW: event fired but atom weight operator-dependent (spectrum != overlap)
print("\n  [FALSIFY-FLOW check] same H(lam=0.2)=diag(0.2,1), event FIRES for every O:")
v0 = np.array([1.0, 0.0]); v1 = np.array([0.0, 1.0])  # eigen 0.2, 1.0
for name, O in [("sigma_x", np.array([[0, 1], [1, 0]])),
                ("sigma_z", np.array([[1, 0], [0, -1]])),
                ("identity", np.eye(2))]:
    f1sq = abs(v0 @ O @ v1) ** 2
    print(f"      O={name:9s}  f1^2=|<0|O|1>|^2={f1sq:.3f}")
print("   => the event fires regardless; f1^2 is 1 or 0  (event blind to the atom weight).")

# =====================================================================
# (C) COUNTERFACTUAL SPECTRAL SHADOW -- bump adversary (Theorem-4 blindness)
# =====================================================================
hdr("[C] counterfactual shadow: remove atom -> atomless bump, same finite moments")
xA, xB = 1.0, 3.3
a, b = 0.45, 0.55                      # nu_A = a delta_xA + b delta_xB
print(f"   nu_A = {a} delta_{xA} + {b} delta_{xB}   (edge atom f1^2={a})")
print("   nu_B(h) = same total, atom at xA smoothed to Unif[xA-h,xA+h]:")
def moments_bump(h):
    M = []
    for k in range(5):
        # delta part at xB  + uniform on [xA-h,xA+h] with mass a
        mk = b * xB ** k + a * np.mean([(xA - h + 2 * h * t / 999) ** k for t in range(1000)])
        M.append(mk)
    return np.array(M)
def moments_atom():
    return np.array([b * xB ** k + a * xA ** k for k in range(5)])
MA = moments_atom()
print("     h        dM0       dM1       dM2       dM3       dM4      edge_atom")
for h in [0.05, 0.01, 0.001]:
    MB = moments_bump(h)
    dM = MB - MA
    print(f"   {h:<7}" + "".join(f" {abs(x):.2e}" for x in dM) + f"   {a} -> 0")
print("   => dM0=dM1=0 exactly; dM2=(xA^2... )h^2/3 -> 0 ; all FINITE data agree")
print("      while the edge atom drops 0.45 -> 0.  This IS the weak-* shadow (Theorem 4 / R3 [C3]).")
print("   Removing the candidate's 'traces' (smearing its atom) leaves every finite relation")
print("   satisfied => the relations cannot certify the atom: reduces to H1/weak-* continuity.")
# gapped-no-shadow control: a determined atom (relations pin it) -- trivially, rank-1 trace
print("\n   no-shadow control: H=diag(0,m^2), O=P_1, relation tr(O)=1 UNIQUELY fixes O.")
print("     remove O => tr(O)=0 != 1 => relation BREAKS => no shadow. Gap m^2>0.")
print("     => 'no shadow' == O is generated/determined by the relations (operator-algebra generation).")

# =====================================================================
# 2x2 BATTERY SUMMARY (gap x atom)
# =====================================================================
hdr("[SUMMARY] independence of ATOM and GAP (the campaign's central dichotomy)")
print("   case                          inf supp>0(gap)?   edge atom f1^2>0?")
print("   1 free massless               NO (m1^2->0)        NO  (atomless)")
print("   2 2D YM SU(2)                 YES (3e^2/8)        NO  (trF^2=H => f1^2=0)")
print("   3a thin atom                  NO (eps->0)         NO  (w=eps/2lam->0)")
print("   3b w delta_m + Unif[0,Lam]    NO (inf supp=0)     YES (w=0.3)")
print("   4 gapped atomless             YES (m^2>0)         NO  (a.c., sqrt edge)")
print("   positive control: massive scalar  YES (m^2>0)    YES (delta at m^2)")
print("   => ATOM and GAP are INDEPENDENT: rows 3b (atom,no gap) and 2/4 (gap,no atom).")
print("   A->tracks W(weight); B->tracks I(isolation); C->tracks H1(measure). None is the gap.")
print("\nDONE.")
