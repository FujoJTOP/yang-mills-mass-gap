"""
tool_nevanlinna.py -- is Nevanlinna-Pick / Caratheodory interpolation on the Stieltjes
function Pi(z)=int dnu/(s-z) a NON-INTEGRAL FORCING lower bound on f1^2 = -Res_{s=1} Pi ?

Claim: NO.  Point evaluations Pi(z_j) are finite INTEGRAL functionals with a bounded kernel
phi_j(s)=1/(s-z_j).  An ATOMLESS bump near s=1 (total mass = the atom weight) mimics an atom's
Pi(z_j) to O(h)->0 as its width h->0  ==>  f1^2=0 is feasible.  NP gives an UPPER bound only.
So NP on point VALUES is in the OPE-RHO integral class and is acid-dirty.

Run:  py -3.12 -I research/unconditional/w7/tool_nevanlinna.py
"""
import numpy as np
from scipy.optimize import linprog

zpts = np.array([0.5 + 0.3j, 0.6 + 0.9j, 1.3 + 0.4j, 2.0 + 0.7j])

# "truth" continuum: dense grid on [2,10], density 0.1  (mass 0.8)
sc = np.linspace(2.0, 10.0, 4000)
wc = 0.1 * (sc[1] - sc[0]) * np.ones_like(sc)
def ev(ss, ws, zs):
    return np.array([np.sum(ws / (ss - z)) for z in zs])

atom_w = 0.5
Pi_atom = ev(sc, wc, zpts) + atom_w / (1.0 - zpts)          # atom at s=1

print("=== (1) atomless bump (f1^2=0) mimics the atom's Pi(z_j) to O(h) ===")
for h in [1e-2, 1e-3, 1e-4, 1e-5]:
    sb = np.linspace(1.0 - h, 1.0 + h, 4001)
    wb = atom_w * np.ones_like(sb) / len(sb)                # total mass = atom_w, NO delta
    Pi_bump = ev(sc, wc, zpts) + ev(sb, wb, zpts)
    print(f"   h={h:8.0e}   max|Pi_bump - Pi_atom| = {np.max(np.abs(Pi_bump - Pi_atom)):.2e}")
print("   => an ATOMLESS measure matches the atom's point evals as h->0  =>  f1^2=0 feasible")

# (2) LP: max atom weight f at s=1 s.t. the same point evals (grid = continuum + atomless fine grid)
g = np.concatenate([np.linspace(0.02, 0.97, 800),
                    np.linspace(1.03, 1.97, 400),
                    np.linspace(2.0, 10.0, 800)])
A = np.array([[1.0 / (s - z) for s in g] for z in zpts])
Ar = np.vstack([A.real, A.imag]); br = np.concatenate([Pi_atom.real, Pi_atom.imag])
c = 1.0 / (1.0 - zpts); cr = np.concatenate([c.real, c.imag])
Aeq = np.hstack([Ar, cr.reshape(-1, 1)])
cobj = np.zeros(Ar.shape[1] + 1); cobj[-1] = -1.0
r = linprog(cobj, A_eq=Aeq, b_eq=br, bounds=[(0, None)] * (Ar.shape[1] + 1), method="highs")
fmax = r.x[-1] if r.success else float("nan")
print("\n=== (2) LP: max atom weight from the SAME point evals ===")
print(f"   f1^2_max (UPPER bound) = {fmax:.4f}   (true atom weight = 0.50)")
print(f"   NP-feasible f1^2 in [0, {fmax:.4f}]  -> includes 0  => NOT forcing")

print("\nCONCLUSION: NP on point VALUES is an INTEGRAL functional (bounded kernel) ->")
print("OPE-RHO applies; atomless bumps mimic any atom, so f1^2=0 is always feasible.")
print("NP yields an UPPER bound and cannot certify f1^2>0 or f1^2=0 (acid-dirty).")
print("The non-integral variant (interpolate ZEROS of Pi) is the A3/read-not-force class.")
print("Blind spot (non-integral FORCING lower-bound atom detector) NOT filled by this family.")
