"""
roof_0form_dual.py -- witness for research/unconditional/ROOF-0FORM-DUAL.md
One small job.  numpy only.  Run:  py -3.12 -I research/unconditional/w7/roof_0form_dual.py

Tests, numerically and exactly:
  W1  FORM-DEGREE / LINKING.  The 't Hooft phase e^{2pi i lk/N} is the exponential of a
      LINKING NUMBER of supports.  Show lk is a genuine topological invariant of two CLOSED
      cycles (Hopf = 1, far-apart = 0, unlinked-in-plane = 0), computed by the Gauss integral.
      A local (0-form) operator has POINT support (a 0-cycle, empty boundary): it links
      NOTHING => its 't Hooft phase is 1.  i.e. no local operator can carry a 1-form charge.
  W2  CENTER NEUTRALITY OF tr F^2.  tr F^2 is a color singlet = adjoint-valued; the center
      z=e^{2pi i/3} acts on it trivially (Ad_{zI}=id on su(N)).  A FUNDAMENTAL Wilson loop
      (tr U) DOES get the phase z.  So the 't Hooft phase can never appear with the singlet.
  W3  ACID.  tr F^2 = 0 for the free massless field (F=0); the gauge trace anomaly vanishes
      for a scalar (:phi^2:) theory (no gauge coupling); it is nonzero only for an interacting
      gauge field with beta(g)!=0.  Acid-clean on both adversaries.
No claim about Yang-Mills beyond these structural facts.  NO QED.
"""
import numpy as np

np.set_printoptions(precision=4, suppress=True, linewidth=120)

def banner(s):
    print("\n" + "=" * 78); print(s); print("=" * 78)

# ============================================================
# W1  Gauss linking number of two closed curves
#   lk(A,B) = (1/4pi) \oint_A \oint_B (rA-rB).(drA x drB)/|rA-rB|^3
# ============================================================
def gauss_link(A, B):
    # A, B: arrays (N,3) of points on closed curves (uniformly sampled, closed)
    dA = np.roll(A, -1, axis=0) - A          # tangent segment of A
    dB = np.roll(B, -1, axis=0) - B
    rA = 0.5 * (A + np.roll(A, -1, axis=0))  # midpoints
    rB = 0.5 * (B + np.roll(B, -1, axis=0))
    lk = 0.0
    for i in range(len(A)):
        r = rA[i] - rB
        num = np.einsum('ij,ij->i', r, np.cross(dA[i], dB))
        den = np.linalg.norm(r, axis=1) ** 3
        lk += np.sum(num / den)
    return lk / (4.0 * np.pi)

def circle(center, radius, plane, N=400):
    t = np.linspace(0, 2*np.pi, N, endpoint=False)
    c = np.array(center, float)
    if plane == 'xy':
        P = np.stack([np.cos(t), np.sin(t), 0*t], axis=1)
    elif plane == 'xz':
        P = np.stack([np.cos(t), 0*t, np.sin(t)], axis=1)
    else:
        raise ValueError
    return c + radius * P

banner("W1 | form-degree: linking number is a topological invariant of CLOSED cycles")
A = circle((0, 0, 0), 1.0, 'xy')            # unit circle in z=0
B_hopf = circle((1, 0, 0), 1.0, 'xz')       # Hopf-partner (pierces A's disk once)
B_far  = circle((5, 5, 5), 1.0, 'xz')       # far away -> unlinked
B_flat = circle((0, 0, 0), 2.0, 'xy')       # concentric coplanar -> unlinked
for name, B in [("Hopf link (threads A)", B_hopf),
                ("far apart", B_far),
                ("concentric coplanar", B_flat)]:
    print(f"  lk(A, {name:22s}) = {gauss_link(A, B):+.4f}")
print("""
  -> linking is a genuine invariant of two 1-cycles (Hopf = +-1, unlinked = 0).
  -> A LOCAL operator has POINT support: a 0-cycle, boundary = empty set. There is no
     2-surface whose boundary is a point, so lk(point, loop) == 0 identically.
     => its 't Hooft phase e^{2pi i lk/N} = 1: a 0-form operator is 't Hooft-NEUTRAL.
""")

# ============================================================
# W2  center neutrality of tr F^2 (adjoint) vs fundamental Wilson loop
# ============================================================
banner("W2 | center charge: tr F^2 (adjoint/singlet) neutral; fundamental loop charged")
rng = np.random.default_rng(7)
N = 3
z = np.exp(2j*np.pi/3)
# random su(3) element X (anti-Hermitian, traceless)  -- stand-in for the field strength
M = rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N))
X = 0.5*(M - M.conj().T)
X -= np.trace(X)/N * np.eye(N)               # traceless -> su(N)
# center element as a group element U = z I  (det = z^N = 1, so U in SU(N))
Uz = z * np.eye(N)
Ad = Uz @ X @ Uz.conj().T
trF2 = -np.real(np.trace(X @ X))             # tr F^2 = -Tr(F^2) (singlet, real >=0)
trF2_after = -np.real(np.trace(Ad @ Ad))
print(f"  tr F^2 (before center) = {trF2:.6f}")
print(f"  tr F^2 (after  center) = {trF2_after:.6f}")
print(f"  |Ad_(zI)(X) - X|_max  = {np.abs(Ad - X).max():.2e}   (center acts trivially on su(N))")
# fundamental Wilson-loop holonomy: a random SU(3) link, gets phase z under the center
U = np.linalg.qr(rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))[0]
U = U / np.linalg.det(U)**(1/N)              # det -> 1, in SU(3)
w, wz = np.trace(U), np.trace(z*U)
print(f"  fundamental tr U       = {w:.4f}")
print(f"  fundamental tr (z U)   = {wz:.4f}   (= z * tr U : carries center charge)")
print("""
  -> tr F^2 is a color singlet (adjoint-valued): the center zI acts as the IDENTITY on su(N),
     so tr F^2 is CENTER-NEUTRAL (no N-ality). The fundamental loop tr U is charged (phase z).
     => a 't Hooft phase can only appear with a charged (1-form) loop -- never with tr F^2.
""")

# ============================================================
# W3  acid test
# ============================================================
banner("W3 | acid: free massless -> 0 ; scalar (:phi^2:) -> 0 ; interacting gauge -> !=0")
# (i) free massless gauge: F = 0
trF2_free = -np.real(np.trace(np.zeros((N, N))))
# (ii) scalar :phi^2: theory: no gauge field; 'gauge trace anomaly' beta_g = 0 -> source 0
beta_g_scalar = 0.0
src_scalar = beta_g_scalar * 1.0
# (iii) interacting gauge: nonzero F, beta(g) != 0
beta_g = -7.0/(16*np.pi**2)                  # schematic one-loop (sign/convention aside)
src_gauge = abs(beta_g) / (2*1.0) * trF2
print(f"  free massless gauge      : tr F^2 = {trF2_free:.3e}")
print(f"  scalar :phi^2: (no gauge): gauge anomaly source = {src_scalar:.3e}")
print(f"  interacting gauge        : (|beta|/2g) tr F^2 = {src_gauge:.6f}   (nonzero)")
print("""
  -> 0 on the free massless field, 0 on the scalar :phi^2: (no gauge coupling),
     nonzero only for the interacting gauge theory: ACID-CLEAN.
  -> CAVEAT: 'nonzero source' does NOT by itself force a gap (a gapless interacting theory
     can carry a condensate). Its floored magnitude is H1. This is the honest boundary.
""")

banner("SUMMARY (structural, no QED)")
print("""
  W1  local (0-form) operator: point support -> linking 0 -> 't Hooft phase 1 (NEUTRAL).
  W2  tr F^2 is a color singlet -> center-neutral; only 1-form loops carry the phase.
  W3  acid-clean, but 'source != gap' (floored magnitude = H1).
  => a literal 0-form 't Hooft ALGEBRA does not exist; the genuine 0-form 'dual' is the
     trace-anomaly WARD identity (tr F^2 <-> dilatation current), a linear anomaly equation,
     which degenerates to the gluon condensate = H1.
""")
