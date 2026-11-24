"""
blast_ym_nonlinear.py  --  witness for BLAST-YM-NONLINEAR.md
The nonlinear-YM attack on H1 via Brascamp-Lieb (positive fiber) + exact rank-one Savvidy.
Run: py -3.12 -I research/unconditional/w7/blast_ym_nonlinear.py
"""
import numpy as np
from scipy import integrate

def banner(s):
    print("\n" + "=" * 72); print(s); print("=" * 72)

# ---------------------------------------------------------------------------
# (1) Brascamp-Lieb on the POSITIVE (log-concave) fiber: Var <= 1/Hess.
#     For a Gaussian (log-concave) fiber with curvature k>0, BL is exact:
#     Var(phi) = 1/k = E[(Hess)^{-1}].  This is the "D_+ <= 0" side.
# ---------------------------------------------------------------------------
banner("(1) BL on the positive fiber: Var <= 1/k  (= 'D_+ <= 0')")
for k in (0.5, 1.0, 3.0):
    # density ~ exp(-k x^2/2), Var = 1/k
    x, w = np.polynomial.hermite_e.hermegauss(80)
    x = x / np.sqrt(k); w = w / w.sum()   # standard normal scaled to Var=1/k
    var = np.sum(w * (x - np.sum(w*x))**2)
    print(f"  k={k:4.1f}: Var={var:.6f}   1/k={1/k:.6f}   Var*(k)={var*k:.6f} (==1, BL tight)")

# ---------------------------------------------------------------------------
# (2) Exact rank-one Savvidy softening: Var(phi^2) for V=-(c/2)phi^2+(lam/4)phi^4
#     saturates at 2/lam, independent of c=gB (the true (condensed) vacuum).
# ---------------------------------------------------------------------------
banner("(2) Exact rank-one Savvidy: Var(phi^2) -> 2/lam, c-independent")
lam = 1.0
def var_phi2(c, lam, n=200001, L=12.0):
    phi = np.linspace(-L, L, n)
    V = -(c/2)*phi**2 + (lam/4)*phi**4
    e = np.exp(-(V - V.min()))
    Z = np.trapezoid(e, phi)
    m1 = np.trapezoid(e*phi**2, phi)/Z
    m2 = np.trapezoid(e*phi**4, phi)/Z
    return m2 - m1**2
print(f"  target 2/lam = {2/lam:.4f}")
for c in (2.0, 5.0, 10.0, 20.0, 50.0):
    v = var_phi2(c, lam)
    print(f"  c=gB={c:5.1f}: Var(phi^2)={v:.4f}   (-> {2/lam:.4f}, bounded, NOT growing with c)")

# ---------------------------------------------------------------------------
# (3) The marginal inequality: G'' = S_bare - D_Sav,  D_Sav = g^2 Var(phi^2) <= 2 g^2/lam.
#     H1 (G''>=0) <=> S_bare > D_Sav.  Show it is a single O(1) competition, near-critical.
# ---------------------------------------------------------------------------
banner("(3) H1 <=> S_bare > D_Sav = g^2 Var(phi^2);  marginal O(1)")
c = 10.0; g = 1.0
vphi2 = var_phi2(c, lam)
Dsav = g*g*vphi2
print(f"  g={g}, Var(phi^2)={vphi2:.4f}  =>  D_Sav = g^2 Var = {Dsav:.4f}   (bound 2g^2/lam = {2*g*g/lam:.4f})")
print("  S_bare  ->  G''=S_bare-D_Sav   (H1 needs G''>=0, i.e. S_bare>D_Sav):")
for Sb in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0):
    G2 = Sb - Dsav
    print(f"    S_bare={Sb:4.1f}: G''={G2:+.4f}   H1_holds={G2>=0}   (margin {G2:+.4f})")
print("  => the margin is O(1) and sign depends on S_bare vs D_Sav: a single marginal inequality.")

# ---------------------------------------------------------------------------
# (4) Acid test: free MASSIVE  -> G''=m^2>0 ; free MASSLESS -> G''=0.
# ---------------------------------------------------------------------------
banner("(4) Acid: free field")
print("  free massless: D_Sav=0, S_bare=0  => G''=0   (H1 fails, as required: no gap for free massless)")
print("  free massive m^2: S_bare=m^2, D_Sav=0 => G''=m^2>0 (gap = m^2)")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
banner("SUMMARY")
print("  BL (positive fiber) gives D_+ <= 0;  exact rank-one Savvidy gives D_Sav=g^2 Var(phi^2)<=2g^2/lam (bounded).")
print("  H1 (G''>=0) reduces to the SINGLE marginal inequality  S_bare > D_Sav = g^2 Var(phi^2).")
print("  The O(1) constant (S_bare vs 2g^2/lam) is the condensate/lam = H1.  Blast does NOT close H1.")
