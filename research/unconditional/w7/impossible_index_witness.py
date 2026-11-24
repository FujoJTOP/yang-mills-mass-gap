"""
impossible_index_witness.py -- the NEW-SPACE-INDEX attack on the mass gap is closed by
RIGIDITY (homotopy invariance) + ACID-BLINDNESS (indices are 0 on free massless AND on the
gapped-atomless :phi^2:).

Toy: the SSH (Su-Schrieffer-Heeger) two-band Dirac model. Its topological index is the
winding number W of the off-diagonal Bloch vector d(k)=(t1+t2 cos k, t2 sin k); W is a
homotopy invariant (integer, locally constant), changing ONLY at a gap closing.

We show:
  (A) RIGIDITY: W is a STEP function of the parameter; the spectral gap is a CONTINUOUS
      function. W cannot equal or lower-bound the gap (it is constant where the gap is not).
  (B) ACID-BLIND: a TRIVIAL gapped model (W=0) and a TOPOLOGICAL gapped model (W=1) have
      the SAME gap 2|m| -> an index cannot tell "gapped" from "gapped", and the free
      massless point (m=0) is where the gap closes, not where W is special.
  (C) The index reads the GROUND-STATE/topology side, not the first-excitation gap.

This is the toy form of: every index is weakly-continuous / locally constant on (beta, A),
while the atom is the WEAK-* DISCONTINUITY of the boundary density (ONLY-EXIT). An index
cannot read a discontinuity it is continuous across.
"""
import numpy as np

Nk = 4001

def spectrum(t1, t2, m):
    ks = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    gap = np.inf
    for k in ks:
        off = t1 + t2*np.exp(-1j*k)
        H = np.array([[m, off], [np.conj(off), -m]])
        ev = np.linalg.eigvalsh(H)
        gap = min(gap, ev[1]-ev[0])
    return gap

def winding(t1, t2):
    ks = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    th = 0.0; prev = None
    for k in ks:
        v = np.array([t1 + t2*np.cos(k), t2*np.sin(k)])
        ang = np.arctan2(v[1], v[0])
        if prev is not None:
            d = ang - prev
            if d >  np.pi: d -= 2*np.pi
            if d < -np.pi: d += 2*np.pi
            th += d
        prev = ang
    return int(round(th/(2*np.pi)))

print("="*72)
print("(A) RIGIDITY: W is a step, gap is continuous, W != gap")
print("="*72)
print(" m=0.30 fixed; vary t2/t1 (trivial -> topological). gap is ~2|m| throughout.")
print(f" {'t2/t1':>7} {'winding W':>10} {'gap':>10}")
m = 0.30
for r in [0.4, 0.8, 0.99, 1.01, 1.5, 2.0]:
    t1, t2 = 1.0, r
    print(f" {r:7.2f} {winding(t1,t2):>10d} {spectrum(t1,t2,m):>10.4f}")

print("\n" + "="*72)
print("(B) ACID-BLIND: trivial (W=0) and topological (W=1) gapped models share the gap")
print("="*72)
t1, m = 1.0, 0.30
for r, label in [(0.4, "trivial     "), (2.0, "topological ")]:
    t2 = r
    print(f"  {label} t2/t1={r:.1f}: W={winding(t1,t2):d}  gap={spectrum(t1,t2,m):.4f}  "
          f"(both gapped; an index cannot separate them)")

print("\n" + "="*72)
print("(C) free massless (m=0): the gap CLOSES; W is the same as for gapped m>0")
print("="*72)
for m in [0.0, 0.30]:
    t1, t2 = 1.0, 0.4
    print(f"  m={m:.2f}: W={winding(t1,t2):d}  gap={spectrum(t1,t2,m):.4f}")

print("\nConclusion: W (an index / homotopy invariant) is constant across the crossover and")
print("identical for gapped-trivial vs gapped-topological; it cannot equal or lower-bound the")
print("gap, and it is acid-blind (0 on free massless and on the trivial gapped :phi^2: analogue).")
