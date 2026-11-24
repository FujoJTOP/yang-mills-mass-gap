"""
tool_potential.py -- witness for TOOL-POTENTIAL.
Tests whether potential theory (log potential U, capacity/equilibrium) can FORCE a LOWER bound
on the atom mass f1^2, or only READ the A3/Stieltjes residue, or is blind to a thin atom.

Concrete measure: nu = f2 * delta_{s1} + (1-f2) * (uniform on [s1,S]).
U(z) = f2 log|z-s1| + (1-f2) * (1/(S-s1)) int_{s1}^{S} log|z-s| ds.
Atom term is the log singularity with coefficient f2 (= the residue, read-not-force).
Continuum term is C_cont(z) = (1-f2)*mean_{[s1,S]} log|z-s| , finite as z->s1.

Three checks:
 (A) READING works: extract f2 from U(s1+i e1)-U(s1+i e2) / log(e1/e2) -> f2 as e->0.
 (B) BLINDNESS: the continuum error floor is O(e/|log(e1/e2)|); to resolve a thin atom f2 the
     resolution must satisfy e << f2, i.e. data precision scales with the UNKNOWN f2 => no uniform
     lower bound. Show error/f2 blows up as f2 shrinks at fixed e.
 (C) CAPACITY / equilibrium: log capacity of [s1,S] is (S-s1)/4 and the equilibrium (arcsine)
     measure depend ONLY on the support => INDEPENDENT of f2 => blind to the atom mass.
 (D) ACID: f2=0 -> no log term -> extracted ~0 (does not misfire on the atomless measure), but a
     THIN atom is also read as ~0 (fails the discriminating requirement "forces f2>0").
Run: py -3.12 -I research/unconditional/w7/tool_potential.py
"""
import numpy as np

def U(z, f2, s1, S, n=400001):
    if S <= s1:
        return f2 * np.log(abs(z - s1))
    s = np.linspace(s1, S, n)
    cont = (1.0 - f2) * np.trapezoid(np.log(np.abs(z - s)), s) / (S - s1)
    atom = f2 * np.log(abs(z - s1)) if f2 > 0 else 0.0
    return atom + cont

def extract(f2, s1, S, e1, e2):
    z1, z2 = s1 + 1j*e1, s1 + 1j*e2
    return (U(z1, f2, s1, S) - U(z2, f2, s1, S)) / np.log(e1/e2)

print("="*70)
print("(A) READING: extract f2 from two complex evaluations, e -> 0")
print("="*70)
s1, S = 1.0, 4.0
for f2 in [0.5, 0.1, 1e-2, 1e-3]:
    est = extract(f2, s1, S, 1e-7, 1e-8)
    print(f"  f2_true={f2:<8.1e}  f2_extract={est:.6f}  err={abs(est-f2):.2e}")
print("  -> converges to f2 (this IS the Stieltjes residue = A3, read-not-force).\n")

print("="*70)
print("(B) BLINDNESS: at FIXED (coarse) resolution e, a thin atom is swamped")
print("="*70)
e1, e2 = 1e-6, 1e-7
for f2 in [0.1, 1e-2, 1e-3, 1e-4]:
    est = extract(f2, s1, S, e1, e2)
    err = abs(est - f2)
    print(f"  [fine e~1e-6] f2={f2:<8.1e}  extract={est:.6f}  err={err:.2e}  err/f2={err/f2:8.2f}")
e1c, e2c = 1e-2, 1e-3
for f2 in [0.1, 1e-2, 1e-3, 1e-4]:
    est = extract(f2, s1, S, e1c, e2c)
    err = abs(est - f2)
    print(f"  [coarse e~1e-2] f2={f2:<8.1e}  extract={est:.6f}  err={err:.2e}  err/f2={err/f2:8.2f}")
print("  -> reading is accurate only at e << f2; at coarse e the continuum error (~e) swamps")
print("     a thin atom; required e << f2 with f2 UNKNOWN => NO uniform lower bound.")
print("  -> AND complex-z U is the A3 object (needs analytic continuation, not measurable).\n")

print("="*70)
print("(C) CAPACITY / equilibrium measure are support-determined => blind to f2")
print("="*70)
cap = (S - s1) / 4.0     # logarithmic capacity of the interval [s1,S]
print(f"  log capacity([{s1},{S}]) = (S-s1)/4 = {cap:.4f}   (independent of f2)")
print("  equilibrium measure on [s1,S] = arcsine density      (independent of f2)")
print("  -> a thin atom at the floor does NOT change the support => capacity BLIND.\n")

print("="*70)
print("(D) ACID: f2=0 (atomless continuum) must give extraction ~0")
print("="*70)
est0 = extract(0.0, s1, S, e1, e2)
print(f"  f2_true=0  f2_extract={est0:.2e}  (atomless measure -> no log term)")
print("  -> acid passes on the atomless continuum, BUT a THIN atom is also read ~0;")
print("     a functional that reads BOTH as ~0 cannot FORCE f2>0. [NG]\n")

print("="*70)
print("VERDICT: potential theory is the Stieltjes/A3 object in disguise.")
print("  Reading = A3 residue (read-not-force); capacity/equilibrium = support-only (blind);")
print("  complex-point evaluations are INTEGRALS (continuous kernel) => OPE-RHO class =>")
print("  no forcing lower bound on a thin atom. [NG]; real yield: a clean WHY.")
print("="*70)
