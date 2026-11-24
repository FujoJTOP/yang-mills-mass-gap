"""
kappa_largen.py -- witness for KAPPA-LARGEN.md  (2026-11-24)

Pure algebra on the rigid-string identity kappa = (m/sqrt(sigma))^2/(8 pi^2).
No Yang-Mills calculation. Shows:
  1. the identity and the band kappa in [0.138,0.193] <-> m/sqrt(sigma) in [3.30,3.90];
  2. the large-N scaling m(N)=m_inf+c/N^2 (1/N^2 for unitary gauge groups),
     kappa_inf = (m_inf/sqrt(sigma))^2/(8 pi^2) ~ band BOTTOM;
  3. "compute kappa at large N" == "compute the large-N gap" (derivative != 0),
     and the flux-tube width ~1/sqrt(sigma) is N-independent => no thin-tube route.
Run:  py -3.12 -I kappa_largen.py
"""
import math

EIGHT_PI2 = 8.0 * math.pi ** 2

def kappa_of_ratio(r):
    return r * r / EIGHT_PI2

def ratio_of_kappa(k):
    return math.sqrt(EIGHT_PI2 * k)

print("=" * 72)
print("1. Rigid-string identity  kappa = (m/sqrt(sigma))^2 / (8 pi^2)")
print("=" * 72)
print(f"   8 pi^2 = {EIGHT_PI2:.4f}")
for r in [3.30, 3.33, 3.50, 3.55, 3.79, 3.90]:
    print(f"   m/sqrt(sigma) = {r:.2f}   ->   kappa = {kappa_of_ratio(r):.4f}")
print()
print("   band kappa in [0.138,0.193]  <=>  m/sqrt(sigma) in [3.30,3.90]")
print(f"   kappa(3.30) = {kappa_of_ratio(3.30):.4f}  (band bottom)")
print(f"   kappa(3.90) = {kappa_of_ratio(3.90):.4f}  (band top)")

print()
print("=" * 72)
print("2. Large-N scaling  m(N) = m_inf + c/N^2   (1/N^2: orientable/planar)")
print("=" * 72)
# representative inherited data [V-num] from CONSTRUCT-LARGEN-MASTER witness
data = [(2, 3.79), (8, 3.33)]
(N1, r1), (N2, r2) = data
A1, A2 = 1.0 / N1 ** 2, 1.0 / N2 ** 2
c = (r1 - r2) / (A1 - A2)
m_inf = r1 - c * A1
print(f"   data [V-num, inherited]: N=2 -> {r1}, N=8 -> {r2}")
print(f"   fit:  m_inf = {m_inf:.4f},  c = {c:.4f}")
for N in [2, 3, 8, 10, 1000]:
    m = m_inf + c / N ** 2
    print(f"   N={N:4d}:  m/sqrt(sigma) = {m:.4f}   kappa = {kappa_of_ratio(m):.4f}")
print(f"   N=inf :  kappa_inf = (m_inf={m_inf:.4f})^2/(8pi^2) = "
      f"{kappa_of_ratio(m_inf):.4f}   <-- band bottom")
print()
print("   => kappa N-independent at leading order; kappa_inf = band bottom;")
print("      1/N^2 corrections raise kappa at finite N toward the SU(3) lattice end.")

print()
print("=" * 72)
print("3. Why large-N computes NOTHING about the VALUE (kappa = gap in disguise)")
print("=" * 72)
print("   d kappa / d(m/sqrt(sigma)) = (m/sqrt(sigma))/(4 pi^2) != 0")
print("   => knowing kappa(N) for all N  <=>  knowing m/sqrt(sigma)(N) for all N;")
print("      large-N kappa_inf = (m_inf/sqrt(sigma))^2/(8 pi^2) = the large-N GAP.")
print()
print("   flux-tube width ~ 1/sqrt(sigma); under 't Hooft scaling sigma ~ Lambda^2")
print("   is N-INDEPENDENT  =>  width is N-independent (NOT thin at large N):")
print("   there is no 1/(sigma w0^2) small parameter, so no thin-tube expansion.")
print()
print("   pure Nambu-Goto (no rigidity) gives kappa = 0 (closed string shrinks,")
print("   unstable) -- the large-N FREE string has no stable kappa>0 closed state.")
print()
print("   Acid test: the whole construction needs a confining string (sigma>0);")
print("   the free massless field and :phi^2: have none => absent (vacuously clean).")
