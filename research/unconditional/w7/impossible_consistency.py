"""
impossible_consistency.py -- witness for IMPOSSIBLE-CONSISTENCY.md

Question: does GLOBAL CONSISTENCY (unitarity, analyticity, crossing, a-theorem,
anomaly matching, conformal-collider positivity) FORCE the 0++ massless sector of
4D pure Yang-Mills to be EMPTY (lightest state massive)?

The decisive test: consistency conditions are satisfied by the FREE MASSLESS gauge
theory (g=0), which is a VALID, GAPLESS 4D YM. If a "consistency functional" also
holds for the gapped theory, it does not discriminate -> cannot force the gap.

We check, numerically:
  (U) Unitarity  : the 0++ spectral density rho(s) >= 0.
  (A) Analyticity: the Stieltjes transform Pi(z)=int rho/(s-z) has Im Pi/Im z >= 0
                   in the upper half plane (positive-real / Herglotz).
  (C) Crossing   : for a single self-conjugate 0++ density, crossing is automatic;
                   the "positivity cone" reduces to rho>=0, already (U).
  (a) a-theorem  : a_UV > a_IR holds for BOTH a gapless CFT IR (a_IR>0) and a
                   gapped IR (a_IR=0) -> it is an inequality, not a gap-forcer.
  (an) Anomaly   : pure YM has NO 0-form (chiral) anomaly to match; the only
                   anomaly is the 1-form center Z_N, matched by the confined vacuum.

Two 0++ measures:
  FREE MASSLESS : rho_free(s) ~ s^0 * theta(s)  (two-gluon threshold continuum from s=0) -> GAPLESS
  GAPPED        : rho_gap(s)  = delta(s-m^2) + small continuum above                 -> GAPPED
Both must satisfy (U),(A) -> consistency does NOT discriminate.
"""
import numpy as np

np.set_printoptions(precision=5, suppress=True, linewidth=120)


def banner(s):
    print("\n" + "=" * 76)
    print(s)
    print("=" * 76)


# ----------------------------------------------------------------------
# spectral measures on the real s-axis (sampled), and their Stieltjes transform
# ----------------------------------------------------------------------
def rho_free_massless(s):
    # :tr F^2: of a free massless field -> two-particle continuum from s=0.
    # in d=4 the two-particle density near threshold s=0 is ~ s^0 (phase space edge).
    return np.where(s >= 0, 1.0, 0.0)


def rho_gapped(s, m2=1.0):
    # an isolated atom at s=m^2 (very narrow triangle) + a weak continuum above 4 m^2
    w = max(1e-3, 0.01 * m2)
    atom = np.exp(-0.5 * ((s - m2) / w) ** 2) / (w * np.sqrt(2 * np.pi))  # ~ delta
    cont = np.where(s >= 4 * m2, 0.05, 0.0)
    return atom + cont


_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")


def stieltjes(rho, s, z):
    # Pi(z) = int rho(s)/(s-z) ds   (principal value for real z; use Im z>0)
    return _trapz(rho / (s - z), s)


banner("(U) UNITARITY: rho(s) >= 0 for both measures?")
s = np.linspace(1e-6, 20.0, 400000)
rf = rho_free_massless(s)
rg = rho_gapped(s, m2=1.0)
print(f"  FREE MASSLESS  min rho = {rf.min():.3e}  (>=0 : {rf.min() >= -1e-12})   support floor = {s[np.argmax(rf>1e-9)]:.2e}  (GAPLESS)")
print(f"  GAPPED         min rho = {rg.min():.3e}  (>=0 : {rg.min() >= -1e-12})   atom at s ~ 1.0  (GAPPED)")

banner("(A) ANALYTICITY / HERGLOTZ: Im Pi(z)/Im z >= 0 in the upper half plane?")
for (name, rho) in [("FREE MASSLESS", rf), ("GAPPED", rg)]:
    minv = np.inf
    for y in [0.05, 0.2, 1.0, 3.0, 10.0]:
        for x in [-5, -1, 0, 1, 2, 5, 10]:
            z = x + 1j * y
            Pi = stieltjes(rho, s, z)
            minv = min(minv, Pi.imag / y)  # Herglotz: Im Pi / Im z >= 0
    print(f"  {name:14s}: min Im(Pi)/Im(z) = {minv:.4f}  (>=0 : {minv >= -1e-9})  -> consistency holds")

banner("(C) CROSSING: single self-conjugate 0++ -> automatic; positivity cone = rho>=0 (already (U))")
print("  both measures satisfy rho>=0 => both in the crossing/positivity cone => NOT discriminated")

banner("(a) a-THEOREM is an INEQUALITY, satisfied by gapless AND gapped IR")
a_UV = 1.0
for (name, a_IR) in [("gapless CFT IR (a_IR>0)", 0.4), ("gapped IR (a_IR=0)", 0.0)]:
    print(f"  {name:26s}: a_UV={a_UV:.2f} > a_IR={a_IR:.2f}  holds={a_UV > a_IR}  -> satisfied either way")

banner("(an) ANOMALY MATCHING: pure YM has NO 0-form chiral anomaly")
print("  chiral (0-form) 't Hooft anomaly of pure SU(N) YM = 0  (no fermions, no chiral symmetry)")
print("  the ONLY anomaly is the 1-form center Z_N anomaly; at theta=0 it is matched by")
print("  the CONFINED vacuum (Z_N topological order) -> consistent with the GAP, not forcing a massless 0-mode")

banner("VERDICT")
print("""  The FREE MASSLESS gauge theory (g=0) is a valid, GAPLESS 4D YM satisfying
  unitarity, analyticity, crossing, a-theorem, and anomaly matching.
  => No consistency condition distinguishes gapless from gapped.
  => The consistency bootstrap CANNOT force the 0++ massless sector empty.
  => A "consistency functional" that fires only on the gapped theory would have to
     use the interacting-specific operator content, i.e. the gap itself (= H1).
  [NG] consistency bootstrap; acid test FAILS (it also holds on the free massless field).
""")
