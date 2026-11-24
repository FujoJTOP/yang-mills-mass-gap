"""
impossible_entanglement.py -- witness for IMPOSSIBLE-ENTANGLEMENT.md (2026-11-24)
Exact free-fermion entanglement spectrum (correlation-matrix technique, numpy only).

THE decisive obstruction the entanglement-gap base must confront:
  the entanglement spectrum tracks the EDGE/topological spectrum (Li-Haldane), NOT the
  bulk gap.  Two SSH phases with the SAME bulk gap 2|t1-t2| have DIFFERENT entanglement
  gaps: the trivial phase (no edge) is entanglement-gapped; the topological phase (zero
  edge modes) is entanglement-GAPLESS (an eigenvalue of the region covariance sits at 1/2).
  So "entanglement spectrum gapped => 0++ bulk gap" is FALSE (it is the edge, not the bulk).

  And the entanglement spectrum is determined by the vacuum state (the measure mu); reading
  its gap IS reading mu.  So the entanglement route does not avoid nu.

Also: (F1) area law does NOT imply a gap -- the free massless chain is gapless (phys gap->0)
with log entanglement entropy; and (F3) the entanglement gap differs in value from the physical
gap even in a free model.
Run:  py -3.12 research/unconditional/w7/impossible_entanglement.py
"""
import numpy as np

def ssh(N, t1, t2):
    H = np.zeros((N, N))
    for i in range(N - 1):
        t = t1 if i % 2 == 0 else t2
        H[i, i + 1] = -t
        H[i + 1, i] = -t
    ev, U = np.linalg.eigh(H)
    order = np.argsort(ev)
    n = np.zeros(N)
    n[order[: N // 2]] = 1.0
    C = U @ np.diag(n) @ U.conj().T
    evs = np.sort(ev)
    # number of near-zero single-particle modes (edge states)
    nedge = int(np.sum(np.abs(evs) < 1e-6))
    return C, evs, nedge

def ent_gap(C, l):
    lam = np.clip(np.linalg.eigvalsh(C[:l, :l]), 1e-12, 1 - 1e-12)
    xi = np.log((1 - lam) / lam)
    S = -np.sum(lam * np.log(lam) + (1 - lam) * np.log(1 - lam))
    return np.min(np.abs(xi)), S

N = 240
print("=" * 84)
print("THE decisive test: same BULK gap, different ENTANGLEMENT gap (edge vs bulk)")
print("=" * 84)
print(f"{'phase':34s} {'bulk_gap':>9s} {'#edge':>6s} {'l':>5s} {'ent_gap':>9s} {'S':>7s}")
for (t1, t2, name, bulk) in [
        (1.0, 0.5, "trivial  t1=1.0,t2=0.5", 2 * abs(1.0 - 0.5)),
        (0.5, 1.0, "TOPOLOGICAL t1=0.5,t2=1.0", 2 * abs(0.5 - 1.0))]:
    C, evs, nedge = ssh(N, t1, t2)
    for l in (N // 2,):
        g, S = ent_gap(C, l)
        print(f"{name:34s} {bulk:9.3f} {nedge:6d} {l:5d} {g:9.5f} {S:7.4f}")
print("  => SAME bulk gap 1.0; trivial ent_gap>0 (gapped), topological ent_gap~0 (edge).")
print("     The entanglement gap reads the EDGE spectrum, not the bulk 0++ gap.")

print()
print("=" * 84)
print("(F1/F3) massless: gapless physics + LOG entropy (not area law); ent_gap != phys gap")
print("=" * 84)
for (t1, t2, name) in [(1.0, 1.0, "massless t1=t2"),
                       (1.0, 0.5, "gapped t1=1,t2=0.5")]:
    C, evs, nedge = ssh(N, t1, t2)
    for l in (N // 4, N // 2):
        g, S = ent_gap(C, l)
        # local physical gap at the Fermi level (bulk; excludes edge states)
        mid = N // 2
        phys = evs[mid + (nedge // 2)] - evs[mid - 1 - (nedge // 2)]
        print(f"  {name:26s} l={l:3d}  phys_gap~{phys:7.4f}  ent_gap={g:8.5f}  S={S:6.4f}")

print()
print("VERDICT of the witness:")
print("  * two phases, same bulk gap 1.0, opposite entanglement gap => the entanglement")
print("    spectrum encodes the EDGE/topological data, NOT the bulk 0++ gap (Li-Haldane).")
print("  * massless: phys gap ~0 with growing (log) entropy => area law is not implied and")
print("    does not imply a gap (Hastings is the converse).")
print("  * the entanglement spectrum is a functional of the vacuum (measure mu) => reading")
print("    its gap IS reading mu => no avoidance of nu; the base reduces to H1.")
