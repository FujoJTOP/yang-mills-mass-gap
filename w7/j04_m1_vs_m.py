"""J04: does any model give m1 > m strictly?  (O = 1/4 trF^2 family)

Units/windows used (printed in every section):
  [A] free theory, O bilinear  -> KL density = two-particle cut, threshold s=4m^2  => m1 = 2m
  [B] 2D YM, O = 1/4 trF^2 = (1/2) tr E^2  -> Casimir diagonal => annihilates vacuum => nu == 0
  [C] M6 = oscillator(omega) (x) free circle(L), m = min(2 omega, 4 pi^2 / L^2)
  [D] single plaquette -> no transfer matrix => both m1, m undefined

Only numpy + fractions.  py -3.12 -I j04_m1_vs_m.py
"""
import numpy as np
from fractions import Fraction as F


def hdr(t):
    print("\n=== " + t + " ===")


# ---------------------------------------------------------------- [A] free theory
def kl_density(s, m):
    """ 1-loop KL density of the bilinear :phi^2: in 4D (shape only):
        rho(s) ~ sqrt(1 - 4m^2/s) * theta(s - 4m^2).  Support starts at s = 4m^2. """
    s = np.atleast_1d(np.asarray(s, float))
    thr = 4.0 * m * m
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(s > thr, np.sqrt(np.clip(1.0 - thr / s, 0.0, None)), 0.0)
    return r


def threshold_scan(m, lo=0.0, hi=None, n=2_000_001):
    """ lowest s on the scan where the density is > 0 (numerical bottom of support) """
    hi = 4.0 * m * m * 1.5 + 1e-6 if hi is None else hi
    s = np.linspace(lo, hi, n)
    r = kl_density(s, m)
    idx = np.nonzero(r > 0)[0]
    return (s[idx[0]] if idx.size else float("nan")), r.max()


hdr("[A] free theory: O = bilinear composite (:phi^2: / two-photon cut), m1 = threshold")
print("   scan window: s in [0, 1.5*4m^2], 2e6 points; density ~ sqrt(1-4m^2/s)*theta")
for m in (1.0, 0.5, 1e-3, 0.0):
    sthr, rmax = threshold_scan(m)
    s_exact = 4.0 * m * m
    print(f"   m = {m:<8.3g} m1 = sqrt(s_thr) = {np.sqrt(sthr):<12.6g} "
          f"2m = {2*m:<12.6g} s_thr_scan = {sthr:<12.6g} s_thr_exact = {s_exact:<12.6g} "
          f"rho_max = {rmax:.4f}")
print("   => m1 = 2m exactly (support is [4m^2,inf)).  m>0: m1=2m>m.  m=0: m1=0=m.")
print("   => atom weight f1^2 = nu({m1^2}) = 0 (a cut, never an atom): ATOM2.md:42")

# ---------------------------------------------------------------- [B] 2D YM
hdr("[B] 2D YM: O = 1/4 trF^2 = (1/2) tr E^2 ; O proportional to the Casimir")
# character basis chi_j of SU(2), j = 0,1/2,1,3/2,2 ;  E^2 has eigenvalue C2 = j(j+1)
js = [F(0), F(1, 2), F(1), F(3, 2), F(2)]
C2 = [j * (j + 1) for j in js]                      # exact rationals
Omat = [[F(1, 2) * C2[n] if n == k else F(0) for k in range(len(js))] for n in range(len(js))]
offdiag = max(abs(Omat[n][k]) for n in range(len(js)) for k in range(len(js)) if n != k)
vac_row = [Omat[0][k] for k in range(len(js))]      # <j=0|O|k>
print("   j        :", "  ".join(f"{str(j):>5}" for j in js))
print("   C2=j(j+1):", "  ".join(f"{str(c):>5}" for c in C2))
print("   O_{nk} = (1/2)C2(n) delta_nk  (diagonal)  max|offdiag| =", offdiag)
print("   <Omega|O|k> =", [str(x) for x in vac_row], " => all zero")
print("   sum_n |<Omega|O|n>|^2 = 0 (exact)  => nu_{O,c} == 0 (empty measure) => m1 = +inf")
g2 = F(1)
E = [g2 / 2 * c for c in C2]                        # H = (g^2/2) E^2 on the chain
print("   E_j = (g^2/2)j(j+1) at g^2=1 :", [str(e) for e in E])
m2d = E[1] - E[0]
print(f"   m(2D) = E_(1/2) - E_0 = {m2d} = {float(m2d):.6f}")
for N in (2, 3):
    beta = F(2 * N, 1) / g2                         # beta = 2N/g^2  (Wilson beta = 2N/g^2)
    lhs = F(N * N - 1, 1) / (2 * beta)              # brief's m = (N^2-1)/(2 beta)
    c2f = F(N * N - 1, 2 * N)                       # C2(fundamental)
    rhs = g2 / 2 * c2f                              # (g^2/2) C2(fund) = single-link gap
    print(f"   N={N}: (N^2-1)/(2 beta)|_(beta=2N/g^2) = {lhs} = {float(lhs):.6f}"
          f"   == (g^2/2)C2(fund) = {rhs} = {float(rhs):.6f}   [{lhs == rhs}]")

# ---------------------------------------------------------------- [C] M6
hdr("[C] M6: H = (p_a^2+p_b^2)/2 + (1/2)omega^2 a^2, b on circle of length L")
omega = 1.0
print("   m = min(2 omega, 4 pi^2 / L^2)  (brief convention);  omega =", omega)
print("   O = quadratric analogue (:a^2:) -> m1 = 2 omega ;  O = a (ATOM2) -> m1 = omega")
print(f"   {'L':>10} {'m':>16} {'m1/m (bilinear)':>18} {'m1/m (linear)':>16}")
rows = []
for L in (8.0, 64.0, 1024.0, 1e6):
    mfull = min(2 * omega, 4 * np.pi ** 2 / L ** 2)
    r2 = 2 * omega / mfull
    r1 = omega / mfull
    rows.append((L, mfull, r2, r1))
    print(f"   {L:>10.0f} {mfull:>16.6e} {r2:>18.6e} {r1:>16.6e}")
Lc = 2 * np.pi / np.sqrt(2 * omega)
Lc1 = 2 * np.pi / np.sqrt(omega)
print(f"   strictness onset (bilinear): L > 2pi/sqrt(2 omega) = {Lc:.4f}"
      f"  ; (linear): L > 2pi/sqrt(omega) = {Lc1:.4f}")
print("   => for every L above the onset: m1 > m STRICTLY, with m->0 as L->inf.")
print("   mechanism: O acts only on the a-factor; the light state is the b zero mode.")

# ---------------------------------------------------------------- [D] single plaquette
hdr("[D] single plaquette")
print("   no time direction => no transfer matrix T => m := -log lambda_2(T) undefined;")
print("   no T => no KL measure nu_O on a Hilbert space => m1 undefined.")
print("   (structural; ATOM2.md:45,49 ; T03-TOY.md sec.1)")

# ---------------------------------------------------------------- summary
hdr("[E] summary table (m1 vs m)")
print("   free (m=0)            m1 = 0        m = 0        equality (zero information)")
print("   free massive composite m1 = 2m      m = m>0      STRICT >   (parity/Wick: <1part|:phi^2:|0>=0)")
print("   2D YM                 m1 = +inf     m = 3/8      empty channel (nu == 0), not a heavy state")
print("   M6                    m1 = 2w (or w) m -> 0      STRICT >   (product sector)")
print("   single plaquette      undefined     undefined    no T")
print("   => m1 > m exists : M6 [T]+[E], free massive composite [T].  Both = 'O misses the light state'.")
