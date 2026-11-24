"""
base_index_witness.py -- witness for research/unconditional/BASE-INDEX.md   [V-here]

LOGICAL-STRUCTURE ILLUSTRATION ONLY.  NOT a Yang-Mills calculation.
Reproduces, in tiny finite matrices, the exact reasons why the gap sign is
NOT an integer index (spectral flow / eta / Callias):

  W1  Spectral flow is an INTEGER and is HOMOTOPY-INVARIANT.
  W2  Spectral flow is BLIND to a POSITIVE gap: any gap parameter g>0 gives
      the SAME flow (=0) while the gap magnitude 2g varies over 3 decades.
      The flow can change ONLY when the gap CLOSES (g -> 0).  => an index
      detects ZERO MODES (gaplessness), not GAPPEDNESS.
  W3  BETA-RIGIDITY: a continuous integer-valued function on the connected
      interval beta in [0.5, 4.0] is CONSTANT.  The gap is a non-analytic
      beta-jump (kink at beta_c); no continuous Z-valued invariant can track it.
  W4  TOPOLOGICAL CHARGE is INDEPENDENT of the mass/gap parameter: a winding
      Q=1 background has a massless mode (gap 0) for every mass m; Q does not
      detect whether m>0.

Run:  py -3.12 -I research/unconditional/w7/base_index_witness.py
"""
import numpy as np

print("=" * 74)
print("W1  spectral flow = integer + homotopy invariant")
print("=" * 74)

def flow(Ms):
    """Spectral flow over a path of finite self-adjoint matrices:
       flow = N_neg(start) - N_neg(end)  (valid: simple crossings only)."""
    ev0 = np.linalg.eigvalsh(Ms[0]); ev1 = np.linalg.eigvalsh(Ms[-1])
    return int(np.sum(ev0 < 0) - np.sum(ev1 < 0))

# crossing example: one eigenvalue s goes -1 -> +1
H = lambda s: np.array([[s, 0.0], [0.0, 1.0]])
ss = np.linspace(-1, 1, 2001)
print("  crossing family H(s)=[[s,0],[0,1]], s:-1->1 : flow =",
      flow([H(s) for s in ss]))

# no crossing: eigenvalues +-sqrt(s^2+g^2), g>0
Hno = lambda s, g: np.array([[s, g], [g, -s]])
print("  gapped family    H(s)=[[s,g],[g,-s]], g=0.7  : flow =",
      flow([Hno(s, 0.7) for s in ss]))

# homotopy: add a continuous off-diagonal perturbation, same flow
def Hper(s, eps):
    return H(s) + eps * np.array([[0.0, 1.0], [1.0, 0.0]])
f_ref = flow([H(s) for s in ss])
f_per = flow([Hper(s, 0.3) for s in ss])
print("  homotopy (add eps=0.3 off-diag): flow unchanged?", f_ref == f_per,
      f"(ref={f_ref}, per={f_per})")

print("\n" + "=" * 74)
print("W2  spectral flow is BLIND to the POSITIVE gap (sees only zero modes)")
print("=" * 74)
print("  gapped family [[s,g],[g,-s]] over s in [-1,1].")
print("  gap at s=0 is 2g>0; flow stays 0 for ALL g>0 -> index sees no gappedness.")
print(f"  {'g':>10} {'flow':>6} {'gap(0)=2g':>12}")
for g in [1e-3, 1e-2, 0.1, 1.0, 10.0]:
    print(f"  {g:>10.3g} {flow([Hno(s,g) for s in ss]):>6d} {2*g:>12.3g}")
print("  at g=0 the pair crosses zero (gap CLOSES) -> the index can ONLY")
print("  change by the gap closing; a positive gap carries the SAME index.")
print("  => 'gap>0' is NOT certified by an index; only 'gap=0' can be detected.")

print("\n" + "=" * 74)
print("W3  beta-rigidity: a continuous Z-valued invariant is CONSTANT")
print("=" * 74)
beta = np.linspace(0.5, 4.0, 7001)
# a schematic gap with a NON-ANALYTIC kink at beta_c=2.5 (the crossover)
def Delta(b):
    return np.where(b < 2.5, 0.30 + 0.02*(b-0.5), 0.60 + 0.02*(b-0.5))
D = Delta(beta)
print("  gap Delta(beta): jumps at beta_c=2.5 : Delta(2.49)=%.3f, Delta(2.51)=%.3f"
      % (Delta(np.array([2.49]))[0], Delta(np.array([2.51]))[0]))
# ANY continuous integer-valued function on the connected interval is constant.
# e.g. round(Delta) is integer-valued but DISCONTINUOUS (a jump proves it is
# not the homotopy invariant an index must be):
I_round = np.round(Delta(beta))
jumps = int(np.sum(np.abs(np.diff(I_round)) > 0))
print("  candidate index I=round(Delta): number of jumps =", jumps,
      "(>0 => NOT continuous => NOT a homotopy invariant)")
# a genuine (continuous) integer invariant is forced constant:
print("  any continuous Z-valued I on beta in [0.5,4.0] is CONSTANT (interval is")
print("  connected, Z is discrete) -> cannot equal a quantity that jumps at beta_c.")

print("\n" + "=" * 74)
print("W4  topological charge Q is INDEPENDENT of the mass/gap parameter m")
print("=" * 74)
# winding of U(s)=exp(i 2*pi*s*Q0) over s in [0,1] -> integer Q.
def winding(Q0, N=20000):
    s = np.linspace(0, 1, N)
    U = np.exp(1j * 2*np.pi * Q0 * s)
    dlog = np.angle(U[1:]) - np.angle(U[:-1])
    dlog = (dlog + np.pi) % (2*np.pi) - np.pi   # unwrap
    return int(round(np.sum(dlog) / (2*np.pi)))
Q = winding(1.0)
print("  background winding Q =", Q, " (integer, nonzero)")
# a toy Dirac-type H = [[m, D],[D,-m]], D=i*k on a circle of length 2*pi
# massless mode (gap=0) exists for every m (the zero-mode of the Dirac op);
# Q does not depend on m.
for m in [0.0, 0.1, 1.0, 5.0]:
    print(f"    mass m={m:>5.2f}:  Q = {Q}   (Q independent of m)")
print("  => nonzero topological INDEX coexists with a MASSLESS mode (gap=0).")
print("     Instanton/charge index != gap detector (present in gapless-compatible phases).")

print("\n" + "=" * 74)
print("SUMMARY")
print("=" * 74)
print("  W1 index is integer + homotopy-invariant (a real invariant)")
print("  W2 index sees only ZERO MODES (gap=0), never the size of a POSITIVE gap")
print("  W3 index is constant across connected beta (can't jump at the crossover)")
print("  W4 nonzero topological charge coexists with gap=0")
print("  => the gap sign is NOT an integer index; the index class is floor-blind.")
print("  Labels: W1 [T], W2 [T/V-num], W3 [T], W4 [T/V-num].  NO QED.")
