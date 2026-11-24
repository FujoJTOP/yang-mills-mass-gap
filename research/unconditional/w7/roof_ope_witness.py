"""
roof_ope_witness.py -- ROOF-OPE numeric witness (2026-11-24)

Question: can OPE + anomalous Ward identity + dispersion relation + positivity
pin the gluon condensate renormalon-free at the crossover?

This is a SCHEMATIC model of the logic, NOT a Yang-Mills calculation.
It reproduces, transparently, the standard SVZ sum-rule structure:

  Pi(Q^2) = C_0(Q^2)  +  A / Q^4  +  ...        (OPE)
  rho(s)  = -s Pi(s)/pi ...   rho(s) >= 0        (reflection positivity / KL)
  Pi(Q^2) = int rho(s) / (s + Q^2) ds            (dispersion)

The perturbative coefficient C_0 has an asymptotic (divergent) series whose
Borel transform has an IR pole at u = 2 (the dimension-4 gluon condensate).
The Borel sum along the real axis is ambiguous; upper/lower contour
prescriptions differ by 2*pi*i*r*e^{-2t}  ->  a real ambiguity +-pi*r*e^{-2t},
which in the standard one-loop correspondence e^{-2t} <-> (Lam^2/Q^2)^2 gives
  delta C_0 = pi*r*Lam^4/Q^4    -- the SAME order as the condensate term A/Q^4.

Run:  py -3.12 -I roof_ope_witness.py
"""
import numpy as np

np.set_printoptions(precision=4, suppress=True)
r    = 1.0      # residue of the Borel pole at u = 2
Lam  = 1.0      # Lambda (units)
beta0 = 1.0     # toy one-loop coefficient

print("=" * 72)
print("PART A -- u=2 IR renormalon: the Borel-sum ambiguity scales as Lam^4/Q^4")
print("=" * 72)
print("Borel  B(u) = r/(2-u)  (pole at u=2 = dim-4 gluon condensate).")
print("S(t)= int_0^inf e^{-u t}B(u) du ; upper vs lower contour differ by")
print("  2*pi*i*r*e^{-2t}  =>  real ambiguity +/- pi*r*e^{-2t}.")
print("1-loop correspondence e^{-2t} <-> (Lam^2/Q^2)^2  =>  ambiguity = pi*r*Lam^4/Q^4.")
print("Condensate term = A*Lam^4/Q^4.  Ratio = pi*r/A, scale-independent.")
print()
print(f"{'Q^2/Lam^2':>10} {'x=Lam^4/Q^4':>13} {'amb=pi*r*x':>12} {'A*x (A=1)':>11} {'amb/(A x)':>11}")
for n in [2, 3, 5, 10, 30, 100]:
    x = 1.0 / n**2
    amb = np.pi * r * x
    Ax  = 1.0 * x
    print(f"{n:10d} {x:13.3e} {amb:12.3e} {Ax:11.3e} {amb/Ax:11.4f}")
print()
print("=> the renormalon ambiguity is a FIXED fraction (pi*r/A) of the")
print("   condensate at EVERY scale. It never shrinks relative to A.")
print("   => there is no renormalon-free 'value' of the condensate from the OPE.")

print()
print("=" * 72)
print("PART B -- two OPE splits of the SAME physical correlator")
print("=" * 72)
# a manifestly positive spectral density (all weights > 0)  ->  rho >= 0
m2 = np.array([0.4, 1.0, 2.5, 6.0, 15.0, 40.0])
w  = np.array([0.30, 0.40, 0.30, 0.50, 0.40, 0.60])
assert np.all(w >= 0)
M0 = np.sum(w)             # leading 1/Q^2 moment
M1 = -np.sum(w * m2)       # 1/Q^4 coefficient  ( = "A_true" from rho )

def Pi(Q2):
    return np.sum(w / (m2 + Q2))

print(f"rho(s)=sum w_j delta(s-m_j^2),  all w_j>=0  (rho>=0 by construction).")
print(f"True 1/Q^4 coefficient of Pi (from rho): A_true = -sum w m^2 = {M1:.6f} (units Lam^4).")
print()
print("OPE-extracted coefficient under prescription s = +/-1:")
print("  A_s = [limit Q2^2 (Pi - M0/Q2)] + s * pi*r*Lam^4")
print(f"{'Q^2/Lam^2':>10} {'Q2^2(Pi-M0/Q2)':>16} {'A(+)':>10} {'A(-)':>10} {'A(+)-A(-)':>12}")
for n in [10, 20, 50, 200, 1000]:
    Q2 = float(n)
    base = Q2**2 * (Pi(Q2) - M0 / Q2)
    Ap = base + np.pi * r * Lam**4
    Am = base - np.pi * r * Lam**4
    print(f"{n:10d} {base:16.6f} {Ap:10.6f} {Am:10.6f} {Ap-Am:12.6f}")
print()
print(f"=> A(+) - A(-) = 2*pi*r*Lam^4 = {2*np.pi*r*Lam**4:.4f} at EVERY Q^2.")
print("   Both prescriptions reproduce the SAME physical rho (and hence the same Pi).")
print("   => the condensate coefficient is only fixed up to +/- pi*r*Lam^4.")

print()
print("=" * 72)
print("PART C -- positivity does NOT pin it")
print("=" * 72)
print("rho >= 0 is the SAME object for both prescriptions (it is the physical")
print("spectral density, prescription-free). Therefore positivity cannot choose")
print("between s=+ and s=-. Any positivity BOUND on the condensate inherits the")
print("+/- pi*r*Lam^4 ambiguity:")
print()
print("  dispersion:  Pi(Q^2) = int rho(s)/(s+Q^2) ds > 0  (for Q^2>0)")
print("  OPE:         Pi = C0(Q^2) + A/Q^4 + ...")
print("  positivity bound:  A > -[C0(Q^2) + s*dC(Q^2)] * Q^4")
print(f"  -> the bound shifts by  s * pi*r*Lam^4 = +/-{np.pi*r*Lam**4:.4f} between prescriptions.")
print()
print("  The ONLY renormalon-free object is rho(s) itself. Bounding the 0++")
print("  atom / support-floor from rho is exactly H1 (finite-moment data are")
print("  floor-blind: Delta_min(K)=0 for all K, SLOT-BOOTSTRAP).")
print()
print("=" * 72)
print("VERDICT: OPE + Ward + dispersion + positivity does NOT pin the")
print("condensate renormalon-free. [NG]  (u=2 IR renormalon; scheme-dependent")
print("split; positivity gives only ambiguous bounds). Exact step: the u=2")
print("subtraction. Clay 0%; no QED.")
print("=" * 72)
