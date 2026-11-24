"""
falsify_flow.py -- witness for FALSIFY-FLOW.md

Does the SPECTRAL FLOW of a deformed H(lam) FORCE the atom f1^2>0?
Claim under test (the campaign's "only-exit"): forcing the atom needs YM-structure or J->inf=H1.
Route to falsify it: spectral flow (a topological integer count of crossings).

Conclusion shown here: NO. Spectral flow is a property of the SPECTRUM (H(lam)); the atom
f1^2 = |<0|O|1>|^2 is an OPERATOR OVERLAP. They are independent data: the same H(lam) with the
same spectral flow supports BOTH f1^2=1 and f1^2=0. Plus rigidity: flow is an integer, constant
on connected components of lam.

Run: py -3.12 research/unconditional/w7/falsify_flow.py
"""
import numpy as np

print("=" * 70)
print("1. SPECTRAL FLOW: a moving level crosses a reference")
print("=" * 70)
# H(lam) = diag(lam-1, lam+1). The (lam-1) level crosses 0 at lam=1.
def crossings_over(grid, track):
    prev = None; n = 0
    for lam in grid:
        e = track(lam)
        if prev is not None and prev < 0 <= e:
            n += 1
        prev = e
    return n

grid = np.linspace(0.0, 2.0, 2001)
flow = crossings_over(grid, lambda lam: lam - 1.0)
print("  H(lam)=diag(lam-1,lam+1); lower level (lam-1) crosses 0 at lam=1.")
print("  spectral flow (upward crossings of 0, lam:0->2) =", flow, "  [expect 1]")

print()
print("=" * 70)
print("2. SAME H(lam), same flow -> atom f1^2 depends on the operator O")
print("=" * 70)
# At the crossing H(0)=diag(-1,1). Eigenstates |0>=[1,0] (E=-1), |1>=[0,1] (E=+1).
# Atom f1^2 = |<0| O |1>|^2 for the interpolating operator O (= 'tr F^2' analog).
v0 = np.array([1.0, 0.0]); v1 = np.array([0.0, 1.0])
sx = np.array([[0.0, 1.0], [1.0, 0.0]])   # off-diagonal: couples 0<->1
Id = np.array([[1.0, 0.0], [0.0, 1.0]])   # identity: does NOT couple 0<->1
SZ = np.array([[1.0, 0.0], [0.0, -1.0]])  # diagonal: does NOT couple 0<->1
def f1sq(O): return float(abs(v0 @ O @ v1) ** 2)

print("  Same H(lam), spectral flow =", flow, "for every O below:")
print("    O = sigma_x : f1^2 =", f1sq(sx), " [atom present]")
print("    O = identity: f1^2 =", f1sq(Id), " [NO atom, yet flow=1]")
print("    O = sigma_z : f1^2 =", f1sq(SZ), " [NO atom, yet flow=1]")

print()
print("=" * 70)
print("3. RIGIDITY: a beta-deformation with NO crossing -> flow constant")
print("=" * 70)
# H_s = diag(1+s, 2+s): both levels move up together; no crossing for s in [0,1].
ss = np.linspace(0.0, 1.0, 1001)
flow_rigid = crossings_over(ss, lambda s: 1.0 + s)
print("  H_s=diag(1+s,2+s), no crossing; flow over s:0->1 =", flow_rigid, "  [expect 0; constant]")

print()
print("=" * 70)
print("VERDICT")
print("=" * 70)
print("  spectral flow = property of H(lam) (the SPECTRUM);")
print("  f1^2 = |<0|O|1>|^2 = property of the OPERATOR OVERLAP.")
print("  They are INDEPENDENT: same flow (=1) arises with f1^2=1 AND with f1^2=0.")
print("  => a nonzero spectral flow does NOT force f1^2>0 (it is a gap/spectrum object,")
print("     not an atom certifier). Rigidity: flow is topological/constant -> no beta-dependent floor.")
print("  The 'only-exit' claim is NOT falsified: spectral flow is a third reading/index-type")
print("  object, blind to the tr F^2 atom.")
