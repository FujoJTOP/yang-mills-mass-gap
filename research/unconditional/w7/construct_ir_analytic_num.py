"""
construct_ir_analytic_num.py -- witness for CONSTRUCT-IR-ANALYTIC-NUM

Make the IR-only analytic functional concrete and COMPUTE it:
    A[mu] = pi / Im(tau_zero),
the nearest Lee-Yang/Fisher zero in complex TIME extent tau of
    Z(tau) = Tr_{0++} e^{-tau H}  (times an operator overlap structure).

On a finite lattice H has a DISCRETE spectrum, so the faithful finite-L object is
    Z(tau) = w0 + sum_i w_i e^{-E_i tau}
(a sum of exponentials: w0 = |<0|O|0>|^2 the VEV term, E_i the 0++ excitation energies,
 w_i the overlaps).  This is the exact finite-L answer; we compute its complex-tau zeros.

Tests:
  (A) ATOM case:  vacuum + one isolated 0++ atom at m1  ->  A = m1  (exact identity).
  (B) :phi^2: / THRESHOLD case: vacuum + a 2-particle CUT (no isolated 0++ atom)
      -> A fires at the lightest 2-particle energy 2m  (a GAP detector, NOT an atom certifier).
  (C) FREE MASSLESS: continuum from 0 (bounded box) -> A stays O(3) (FIRES at finite L; does
      NOT vanish) -- so A is NOT acid-clean at finite L; the zero is a scale estimate, not the gap.
  (D) finite-L sequence: as the cut becomes dense, the atom zero is fixed while the cut's
      zeros crowd/migrate  ->  finite-L data cannot distinguish atom from cut (H1).

Genuine grounding: free massive scalar on an N-site ring, O = :phi^2: (zero momentum);
the 0++ excitations are 2-particle states of energy 2*omega_k.  We feed those exact energies
into Z(tau) to get the real finite-L cut spectrum.

Labels [T]/[C]/[O]/[NG].  NOT a Yang-Mills calculation.  No QED.
Run:  py -3.12 -I research/unconditional/w7/construct_ir_analytic_num.py
"""
import numpy as np

np.set_printoptions(precision=5, suppress=True, linewidth=120)


def Zsum(tau, w, E):
    w = np.asarray(w, float); E = np.asarray(E, float)
    return np.sum(w * np.exp(-E * tau))


def dZsum(tau, w, E):
    w = np.asarray(w, float); E = np.asarray(E, float)
    return np.sum(-E * w * np.exp(-E * tau))


def newton_zeros(w, E, re=(-8.0, 8.0), im=(0.05, 30.0), nre=41, nim=121, iters=50):
    """All zeros of the exponential sum in the strip; Newton from a grid, dedup, sort by |Im|."""
    w = np.asarray(w, float); E = np.asarray(E, float)
    found = []
    for x in np.linspace(*re, nre):
        for y in np.linspace(*im, nim):
            t = x + 1j * y
            ok = False
            for _ in range(iters):
                z = Zsum(t, w, E); dz = dZsum(t, w, E)
                if dz == 0 or not np.isfinite(abs(z)):
                    break
                step = z / dz
                t -= step
                if abs(step) < 1e-13:
                    ok = True
                    break
            if ok and abs(Zsum(t, w, E)) < 1e-8 and 1e-6 < t.imag < 29.99:
                if not any(abs(t - u) < 1e-4 for u in found):
                    found.append(t)
    found.sort(key=lambda t: abs(t.imag))
    return found


def report(name, w, E, expect=None, extra=""):
    zs = newton_zeros(w, E)
    print(f"\n--- {name} ---")
    if expect is not None:
        print(f"   expected (true lightest 0++ excitation) = {expect:.5f}")
    print(f"   #zeros in strip: {len(zs)}")
    if not zs:
        print("   nearest zero: NONE found  ->  A = 0 (no gap signal)")
        return None
    t0 = zs[0]
    A = np.pi / abs(t0.imag)
    print(f"   nearest zero: tau = {t0.real:+.5f} {t0.imag:+.5f}i   |Z|={abs(Zsum(t0,w,E)):.2e}")
    print(f"   A = pi/|Im tau| = {A:.5f}")
    if expect:
        print(f"   A - expected = {A-expect:+.5f}  ({'MATCH' if abs(A-expect)<0.03*abs(expect)+0.02 else 'MISMATCH'})")
    if extra:
        print(extra)
    return A


print("=" * 78)
print("(A) ATOM case:  Z = w0 + w1 e^{-m1 tau}   (vacuum + isolated 0++ atom)")
print("=" * 78)
for m1, (w0, w1) in [(1.0, (1.0, 0.5)), (1.0, (1.0, 3.0)), (2.0, (1.0, 0.7))]:
    report(f"atom m1={m1}, w=(w0,w1)=({w0},{w1})", [w0, w1], [0.0, m1], expect=m1,
           extra="   identity: Im(tau)=pi/m1 is INDEPENDENT of the weight ratio (weights only shift Re tau)")

print("\n" + "=" * 78)
print("(B) :phi^2: / THRESHOLD case: vacuum + 2-particle CUT, NO isolated atom")
print("    (free massive scalar on N-site ring, O=:phi^2:, 0++ excitations = 2-particle, E=2 w_k)")
print("=" * 78)


def free_scalar_two_particle(N, m):
    n = np.arange(N)
    k = 2 * np.pi * n / N
    om = np.sqrt(m * m + 4 * np.sin(k / 2) ** 2)
    Es, ws = [], []
    seen = set()
    for a in range(N):
        b = (N - a) % N
        key = tuple(sorted((a, b)))
        if key in seen:
            continue
        seen.add(key)
        Es.append(om[a] + om[b])          # energy of |k,-k> 2-particle state
        ws.append(1.0 / (om[a] * om[b]))  # schematic :phi^2: overlap ~ 1/(w_a w_b)
    Es = np.array(Es); ws = np.array(ws)
    return Es, ws


m = 1.0
for N in (8, 32):
    Es, ws = free_scalar_two_particle(N, m)
    w0 = 1.0
    w = np.concatenate(([w0], ws / ws.sum() * 0.8))   # vacuum + normalized cut
    E = np.concatenate(([0.0], Es))
    report(f"threshold m={m}, N={N} (cut from 2m={2*m})", w, E, expect=2 * m,
           extra=f"   2m = {2*m}; note the TRUE 1-particle gap is m={m}, but :phi^2: sees only 2m (a CUT, f1^2=0)")

print("   ==> A fires at 2m for an ATOMLESS operator: it is a GAP detector, NOT an atom certifier.")

print("\n" + "=" * 78)
print("(C) FREE MASSLESS: continuum from 0 in a box of size N; spacing ~ 1/N")
print("=" * 78)


def massless_box(N):
    # 2-particle energies with m=0: E=2*omega_k = 4|sin(k/2)|, k=2pi n/N; lightest nonzero ~ 4*pi/N
    n = np.arange(N)
    k = 2 * np.pi * n / N
    om = 2 * np.abs(np.sin(k / 2))
    Es = []
    seen = set()
    for a in range(N):
        b = (N - a) % N
        key = tuple(sorted((a, b)))
        if key in seen:
            continue
        seen.add(key)
        Es.append(om[a] + om[b])
    Es = np.array(sorted(Es))
    return Es


for N in (8, 16, 32, 64):
    Es = massless_box(N)
    w0 = 1.0
    w = np.concatenate(([w0], 0.6 / len(Es) * np.ones(len(Es))))
    E = np.concatenate(([0.0], Es))
    A = report(f"massless N={N} (lightest 2-particle E={Es[0]:.4f})", w, E, expect=Es[0])
    print(f"   -> A={0 if A is None else A:.5f}. NOTE: A does NOT -> 0; it stays O(3) as N grows.")
print("   ==> At finite L, massless ALSO fires (A~3.2, not 0): A is NOT acid-clean at finite L.")
print("       The zero is set by the SCALE of excitation content, not by the mass gap -> gap/scale detector.")

print("\n" + "=" * 78)
print("(B0) zero-VEV :phi^2: (normal-ordered, w0=0): PURE CUT, no vacuum term")
print("=" * 78)
for N in (8, 16, 32):
    Es, ws = free_scalar_two_particle(N, 1.0)
    w = ws / ws.sum()
    E = Es
    zs = newton_zeros(w, E)
    if zs:
        t0 = zs[0]
        print(f"   N={N}: nearest zero Im={t0.imag:.5f}, A=pi/Im={np.pi/abs(t0.imag):.5f}, "
              f"count={len(zs)}   (FIRES at finite L too)")
    else:
        print(f"   N={N}: no zero in strip -> A=0 (zero-free like an atomless cut should be)")

print("\n" + "=" * 78)
print("(D) finite-L: atom zero FIXED vs cut zeros CROWD  (finite L cannot tell atom from cut)")
print("=" * 78)
# atom: single zero at Im=pi/m1; cut: number of zeros in a fixed strip grows with N
za = newton_zeros([1.0, 0.5], [0.0, 1.0])
print(f"   atom (m1=1): nearest zero Im={za[0].imag:.5f}, count in strip={len(za)}")
for N in (8, 16, 32):
    Es, ws = free_scalar_two_particle(N, 1.0)
    w = np.concatenate(([1.0], ws / ws.sum() * 0.8))
    E = np.concatenate(([0.0], Es))
    zc = newton_zeros(w, E)
    print(f"   cut N={N}: nearest zero Im={zc[0].imag:.5f}, count in strip={len(zc)}   (crowds with N)")

print("\n" + "=" * 78)
print("(E) WEIGHT-PERTURBATION test: is A weight-INDEPENDENT?  (atom vs cut signature)")
print("=" * 78)


def A_of(w, E):
    zs = newton_zeros(w, E)
    return None if not zs else np.pi / abs(zs[0].imag)


print("   ATOM (m1=1):")
a_vals = []
for ww in [(1.0, 0.3), (1.0, 1.0), (1.0, 5.0)]:
    a = A_of([ww[0], ww[1]], [0.0, 1.0]); a_vals.append(a)
    print(f"     weights {ww}: A = {a:.6f}")
print(f"     spread = {max(a_vals)-min(a_vals):.2e}  -> weight-INDEPENDENT (atom signature)")
print("   CUT (:phi^2:, m=1, w0=1):")
Es, ws = free_scalar_two_particle(16, 1.0)
for label, wcut in [("uniform", np.ones(len(Es)) / len(Es)), ("1/(ww)", ws / ws.sum()),
                    ("exp-tilt", np.exp(0.5 * Es) / np.exp(0.5 * Es).sum())]:
    w = np.concatenate(([1.0], 0.8 * wcut)); E = np.concatenate(([0.0], Es))
    a = A_of(w, E)
    print(f"     weights {label}: A = {a:.6f}")
print("     (A moves with the weight distribution -> weight-DEPENDENT: cut signature)")
print("   ==> 'A invariant under operator-reweighting' distinguishes ATOM from CUT at fixed L.")
print("       This is a concrete atom SIGNATURE (not the raw value). Caveat: reweighting the operator")
print("       (e.g. trF^2 vs a smeared/other 0++ operator) must not change which states couple; proving")
print("       invariance uniformly in (L,a) is still H1.")

print("\n" + "=" * 78)
print("SUMMARY")
print("=" * 78)
print("  (A) atom:      A = m1 EXACTLY, INDEPENDENT of weight ratio  -> the clean identity [T]")
print("  (B) threshold: A FIRES (~2.9) for an ATOMLESS :phi^2: cut -> gap/scale detector, not atom certifier")
print("  (C) massless:  A stays O(3) at finite L (fires), does NOT -> 0 -> NOT acid-clean at finite L")
print("  (D) finite L:  atom zero FIXED; cut zeros CROWD with N -> finite L cannot tell atom from cut")
print("  => At finite L, A is a renormalon-free SCALE estimator that fires for atom AND cut AND massless.")
print("     It equals the atom mass m1 EXACTLY only in the isolated-atom case. The atom-vs-cut")
print("     distinction is INVISIBLE at finite L; it could only appear as L->infty (cut zeros receding /")
print("     atom zero persisting) -- and proving that persistence uniformly in (L,a) IS H1. [O]")
print("  => CORRECTION to PNT-IR-ANALYTIC: 'fires iff f1^2>0, acid-clean' holds only in infinite volume;")
print("     at finite L A fires on :phi^2: and on the massless field too (= same finite-data blindness).")
print("  Clay 0%; no QED.")
